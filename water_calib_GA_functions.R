### WATER CALIBRATION ALGORITHM - FUNCTIONS
### Author: Kevin J. Wolz

read_profiles <- function(path, profile, hisafe = TRUE) {
  profile.files <- list.files(path, pattern = paste0(profile, "\\.txt"), full.names = TRUE, recursive = TRUE)
  out <- purrr::map_df(profile.files, read.table,
                       header           = TRUE,
                       sep              = "\t",
                       stringsAsFactors = FALSE,
                       na.strings       = c("NA", "error!", "NaN", "-9999", "Infinity"),
                       encoding         = "latin1") %>%
    dplyr::as_tibble()

  if(hisafe) {
    out <- out %>%
      dplyr::mutate(id  = as.numeric(str_remove(SimulationName, paste0("gen_", GEN, "-")))) %>%
      dplyr::mutate(gen = GEN) %>%
      dplyr::select(gen, id, dplyr::everything())
  }
  return(out)
}

write_log <- function(x, sep = "\n") cat(x, file = log.file, sep = sep, append = TRUE)

copy_template <- function(id, overwrite = TRUE) {
  template.folder <- paste0(template.path, "hisafe")

  dum <- file.copy(from      = template.folder,
                   to        = gen.folder,
                   recursive = TRUE,
                   overwrite = overwrite)

  dum <- file.rename(paste0(gen.folder, "hisafe"),
                     paste0(gen.folder, "gen_", GEN, "-", id))
  dum <- file.rename(paste0(gen.folder, "gen_", GEN, "-", id, "/hisafe.sim"),
                     paste0(gen.folder, "gen_", GEN, "-", id, "/gen_", GEN, "-", id, ".sim"))
}

move_lai <- function(id) {
  dum <- file.copy(from = paste0(output.path, "lai.obs"),
                   to   = paste0(gen.folder, "gen_", GEN, "-", id))
}

round_params <- function(params) {
  for(i in 1:length(params)){
    params[i] <- round(params[i] / PARAMS$sig.level[i]) * PARAMS$sig.level[i]
  }

  return(params)
}

edit_crop_parameters <- function(id, param) {
  file.path <- paste0(gen.folder, "gen_", GEN, "-", id, "/cropSpecies/", CROP, ".plt")
  crop.file <- sim <- scan(file     = file.path,
                           what     = "character",
                           encoding = "latin1",
                           sep      = "\n",
                           quiet    = TRUE)

  for(i in 1:length(param)) {
    location <- grep(PARAMS$param.name[i], crop.file)
    crop.file[location] <- gsub(" = .*", paste(" =", param[[i]]), crop.file[location])
  }

  write(crop.file, file = file.path)
}

write_script <- function(ids) {

  write_line <- function(x) cat(x, file = batch.script, sep = "\n", append = TRUE)

  batch.script <- file(description = paste0(gen.folder, "gen_", GEN, ".sh"),
                       open        = "wb",
                       encoding    = "UTF-8")
  cat("", file = batch.script, sep = "", append = FALSE)
  write_line("#!/bin/sh")
  write_line(paste0("#SBATCH --array=", min(ids), "-", max(ids), "%28"))
  write_line("#SBATCH --account=hisafe")
  write_line("#SBATCH --partition=defq")
  #write_line("#SBATCH --account=f_meso")
  #write_line("#SBATCH --partition=fmuse1")
  write_line("#SBATCH --mail-type=NONE")
  write_line("#SBATCH --mail-user=wolzkevin@gmail.com")
  write_line("module purge")
  write_line("module load jre/jre.8_x64")
  write_line("cd /nfs/work/hisafe/Capsis4")
  write_line(paste0("sh capsis.sh -p script safe.pgms.ScriptGen ",
                    gen.folder,
                    "gen_", GEN, "-$SLURM_ARRAY_TASK_ID", "/",
                    "gen_", GEN, "-$SLURM_ARRAY_TASK_ID.sim water_calib_defaults"))
  close(batch.script)
}

check_completion <- function(ids, path) {
  num.completed <- list.files(path       = path,
                              pattern    = "session\\.txt",
                              recursive  = TRUE,
                              full.names = TRUE) %>%
    purrr::map(scan,
               what = "character",
               quiet = TRUE) %>%
    purrr::map_lgl(function(x) any(stringr::str_detect(x, pattern = "Duration"))) %>%
    sum()

  return(num.completed == length(ids))
}

summarize_fitness <- function(fitness) {
  fitness.summary <- fitness %>%
    t() %>%
    dplyr::as_tibble()

  n.sims <- fitness.summary %>%
    summarize(n = n()) %>%
    .$n

  fitness.summary <- fitness.summary %>%
    tidyr::gather(key = "metric", value = "fitness") %>%
    dplyr::group_by(metric) %>%
    dplyr::summarize_all(c("mean", "median", "min", "max", "sd")) %>%
    dplyr::mutate(gen = GEN) %>%
    select(gen, metric, dplyr::everything()) %>%
    dplyr::mutate(n = n.sims)
  return(fitness.summary)
}

crop_water_GA_fitness <- function(population) {
  pop <- do.call(rbind, population) %>%
    as.data.frame()
  names(pop) <- PARAMS$param.name

  ids <- 1:nrow(pop)

  ### DETERMINE IDS ALREADY RUN
  # write_log("----- determine ids already run")
  # out <- pop %>%
  #   dplyr::left_join(pop.database, by = PARAMS$param.name)
  #
  # already.run <-

  pop.to.store <- pop %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(gen = GEN) %>%
    dplyr::mutate(id = ids)

  ### CREATE FOLDER STRUCTURE
  write_log("----- create folder structure")
  dir.create(gen.folder, showWarnings = FALSE)
  purrr::walk(ids, copy_template)
  purrr::walk(ids,  move_lai)

  ### EDIT CROP PARAMETERS
  write_log("----- edit crop parameters")
  purrr::walk2(ids,
               split(as.matrix(pop[ids,]), seq(nrow(pop[ids,]))),
               edit_crop_parameters)

  ### CREATE BATCH FILE
  write_log("----- create batch file")
  write_script(ids)

  ### RUN SUMULATIONS
  write_log("----- run simulations")
  startSims <- proc.time()[3]
  system(paste0("sbatch ", gen.folder, "gen_", GEN, ".sh"))

  ### MONITOR SIMULATIONS
  incomplete <- TRUE
  while(incomplete) {
    Sys.sleep(1)
    complete <- check_completion(ids, gen.folder)
    if(complete) incomplete <- FALSE
  }

  elapsedSims = round(proc.time()[3] - startSims)
  write_log(paste0("----- done with simulations (", elapsedSims, " seconds)"))

  ### READ SIMULATION RESULTS
  write_log("----- read simulations")
  sim.data <- read_profiles(gen.folder, "voxelsOptim")
  HISAFE <- sim.data %>%
    dplyr::filter(Date %in% GROWTH.DATES) %>%
    dplyr::select(gen, id, Date, z, cropWaterUptake) %>%
    dplyr::rename(hisafe = cropWaterUptake)

  ### CALCULATE FITNESS VALUES
  write_log("----- calculate fitness")

  joined.data <- HISAFE %>%
    dplyr::left_join(STICS, by = c("Date", "z"))

  layer.rmse <- joined.data %>%
    dplyr::filter(z <= MAX.DEPTH) %>%
    dplyr::mutate(sqdif = (hisafe - stics) ^ 2) %>%
    dplyr::group_by(gen, id, z) %>%
    dplyr::summarize(rmse = sqrt(mean(sqdif))) %>%
    dplyr::ungroup() %>%
    tidyr::spread(key = "z", value = "rmse")

  total.uptake <- joined.data %>%
    dplyr::mutate(dif = hisafe - stics) %>%
    dplyr::group_by(gen, id) %>%
    dplyr::summarize(dif = abs(sum(dif))) %>%
    dplyr::ungroup()

  dat.vol <- dat <- layer.rmse %>%
    dplyr::left_join(total.uptake, by = c("gen", "id"))

  fitness <- dat %>%
    dplyr::select(-gen, -id) %>%
    as.matrix() %>%
    t()

  dat <- dat %>%
    dplyr::left_join(pop.to.store, by = c("gen", "id"))

  ### WRITE OUT RESULTS
  write_log("----- write out population data")
  readr::write_csv(dat, paste0(output.path, "GA_data_", GEN, ".csv"))

  ## PLOT ALL SIMULATIONS FROM THIS GENERATION
  test.dat <- joined.data %>%
    dplyr::mutate(Date = lubridate::dmy(Date)) %>%
    dplyr::mutate(Month = lubridate::month(Date)) %>%
    dplyr::filter(Month >= 4, Month <= 6) %>%
    dplyr::group_by(gen, id, Date) %>%
    dplyr::summarize(hisafe = sum(hisafe),
                     stics  = sum(stics))

  test.dat$season <- 2004
  test.dat$season[which(test.dat$Date >= lubridate::dmy("1/10/2004"))] <- 2005
  test.dat$season[which(test.dat$Date >= lubridate::dmy("1/10/2005"))] <- 2006

  stics.ref <- test.dat %>%
    dplyr::ungroup() %>%
    dplyr::filter(id == 1) %>%
    dplyr::select(-gen, -id, -hisafe)

  pop.comp.plot <- ggplot(test.dat, aes(x = Date, y = hisafe, color = factor(id), group = id)) +
    facet_wrap(~season, ncol = 1, scales = "free_x") +
    geom_line(size = 0.5) +
    guides(color = FALSE) +
    labs(y = "cropWaterUptake") +
    geom_line(data = stics.ref, aes(x = Date, y = stics),
              inherit.aes = FALSE, color = "black") +
    theme_bw()
  dum <- ggsave_fitmax(paste0(output.path, "Population_Comp_", GEN, ".pdf"), pop.comp.plot)


  ## PLOT BEST SOLUTION FROM THIS GENERATION
  best.id <- dat.vol %>%
    dplyr::select(-gen, -id) %>%
    as.matrix() %>%
    apply(1, prod) %>%
    which.min()

  best.data.all.layers <- joined.data %>%
    dplyr::left_join(OLD.WINNER, by = c("Date", "z")) %>%
    dplyr::filter(id == best.id) %>%
    dplyr::filter(Date %in% GROWTH.DATES) %>%
    tidyr::gather(hisafe, stics, old.winner, key = "sim", value = "waterUptake") %>%
    dplyr::mutate(Date = lubridate::dmy(Date))

  best.data.all.layers$season <- 2004
  best.data.all.layers$season[which(best.data.all.layers$Date >= lubridate::dmy("1/10/2004"))] <- 2005
  best.data.all.layers$season[which(best.data.all.layers$Date >= lubridate::dmy("1/10/2005"))] <- 2006

  best.data <- best.data.all.layers %>%
    dplyr::filter(z <= 0.9)

  best.total <- best.data.all.layers %>%
    dplyr::group_by(gen, id, Date, sim, season) %>%
    dplyr::summarize(waterUptake = sum(waterUptake)) %>%
    dplyr::mutate(z = "total")

  best.data.all <- best.data %>%
    dplyr::mutate(z = as.character(z)) %>%
    dplyr::bind_rows(best.total)

  best.plot <- ggplot(best.data.all, aes(x = Date, y = waterUptake, color = sim)) +
    ggplot2::labs(x = "Date", y = "Water uptake", color = NULL) +
    ggplot2::facet_grid(z~season, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::scale_color_manual(values = c("green", "red", "black"))
  dum <- ggsave_fitmax(paste0(output.path, "Best_Offspring_Comp_", GEN, ".pdf"), best.plot)

  return(fitness)
}

ggsave_fitmax <- function(filename,
                          plot,
                          maxheight = 7,
                          maxwidth  = maxheight,
                          units     = "in", ...) {
  if(is.null(plot)) return(FALSE)
  dims = get_dims(ggobj     = plot,
                  maxheight = maxheight,
                  maxwidth  = maxwidth,
                  units     = units)
  ggplot2::ggsave(filename = filename,
                  plot   = plot,
                  height = dims$height,
                  width  = dims$width,
                  units  = units, ...)
}

get_dims <- function(ggobj,
                     maxheight,
                     maxwidth = maxheight,
                     units    = "in", ...) {

  # Internal helper function:
  # Treat all null units in a unit object as if they were inches.
  # This is a bad idea in gneral, but I use it here as a workaround.
  # Extracting unit names from non-atomic unit objects is a pain,
  # so questions like "which rows of this table layout have null heights?"
  # are hard to answer. To work around it, I exploit an (undocumented!)
  # quirk: When calculating the size of a table layout inside a Grid plot,
  # convertUnit(...) treats null units as zero.
  # Therefore
  #	(convertHeight(grob_height, "in", valueOnly=TRUE)
  #	- convertHeight(null_as_if_inch(grob_height), "in", valueOnly=TRUE))
  # does the inverse of convertUnit: It gives the sum of all *null* heights
  # in the object, treating *fixed* units as zero.
  #
  # Warning: I repeat, this approach ONLY makes any sense if
  #	convertUnit(unit(1, "null"), "in", "x", valueOnly=T) == 0
  # is true. Please check that it is before calling this code.
  .null_as_if_inch = function(u){
    if(!grid::is.unit(u)) return(u)
    if(is.atomic(u)){
      if("null" %in% attr(u, "unit")){
        d = attr(u, "data")
        u = unit(
          x=as.vector(u),
          units=gsub("null", "in", attr(u, "unit")),
          data=d)
      }
      return(u)
    }
    if(inherits(u, "unit.arithmetic")){
      l = .null_as_if_inch(u$arg1)
      r = .null_as_if_inch(u$arg2)
      if(is.null(r)){
        args=list(l)
      }else{
        args=list(l,r)
      }
      return(do.call(u$fname, args))
    }
    if(inherits(u, "unit.list")){
      return(do.call(grid::unit.c, lapply(u, .null_as_if_inch)))
    }
    return(u)
  }

  if(inherits(ggobj, "ggplot") && !isTRUE(ggobj$respect) &&
     is.null(ggobj$theme$aspect.ratio) && is.null(ggobj$coordinates$ratio) &&
     is.null(ggplot2::theme_get()$aspect.ratio)) {
    return(list(height = maxheight, width = maxwidth))
  }

  tmpf = tempfile(pattern = "dispos-a-plot", fileext = ".png")
  png(filename = tmpf,
      height   = maxheight,
      width    = maxwidth,
      units    = units,
      res      = 120, ...)

  on.exit({
    dev.off()
    unlink(tmpf)
  })

  if (inherits(ggobj, "ggplot")) {
    g = ggplot2::ggplotGrob(ggobj)
  } else if (inherits(ggobj, "gtable")) {
    g = ggobj
  } else {
    stop("Don't know how to get sizes for object of class ", deparse(class(ggobj)))
  }

  stopifnot(grid::convertUnit(unit(1, "null"), "in", "x", valueOnly = TRUE) == 0)
  known_ht = sum(grid::convertHeight(g$heights, units, valueOnly = TRUE))
  known_wd = sum(grid::convertWidth(g$widths,   units, valueOnly = TRUE))
  free_ht  = maxheight - known_ht
  free_wd  = maxwidth  - known_wd
  all_null_rowhts = (grid::convertHeight(.null_as_if_inch(g$heights),
                                         "in", valueOnly = TRUE) - grid::convertHeight(g$heights,
                                                                                       "in", valueOnly = TRUE))
  all_null_colwds = (grid::convertWidth(.null_as_if_inch(g$widths),
                                        "in", valueOnly = TRUE) - grid::convertWidth(g$widths,
                                                                                     "in", valueOnly = TRUE))
  null_rowhts = all_null_rowhts[all_null_rowhts > 0]
  null_colwds = all_null_colwds[all_null_colwds > 0]
  panel_asps = matrix(null_rowhts, ncol = 1) %*% matrix(1 / null_colwds, nrow = 1)
  max_rowhts = free_ht/sum(null_rowhts) * null_rowhts
  max_colwds = free_wd/sum(null_colwds) * null_colwds
  rowhts_if_maxwd = max_colwds[1] * panel_asps[, 1]
  colwds_if_maxht = max_rowhts[1]/panel_asps[1, ]
  height = min(maxheight, known_ht + sum(rowhts_if_maxwd))
  width  = min(maxwidth,  known_wd + sum(colwds_if_maxht))
  return(list(height = height, width = width))
}
