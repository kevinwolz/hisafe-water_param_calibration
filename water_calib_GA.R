### WATER CALIBRATION ALGORITHM - MAIN
### Author: Kevin J. Wolz

GA.id <- 1
CROP <- "durum-wheat-restinclieres" #c("durum-wheat-restinclieres", "weed-restinclieres", "rape", "winter-pea")
PATH <- "/lustre/lecomtei/hisafe-water_param_calibration/"
#PATH <- "./"

##### USER PARAMETERS #####
source(paste0(PATH, "water_calib_GA_functions.R"))

input.path    <- paste0(PATH, "input/")
run.path      <- paste0(PATH, "output", GA.id, "/")
output.path   <- paste0(run.path, CROP, "/")
gen.path      <- paste0(output.path, "generations/")
template.path <- paste0(run.path, "templates/", CROP, "/")

dir.create(gen.path, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)
library(purrr)
library(stringr)
library(GGally)
library(ecr)

##### GA CONTROL PARAMETERS #####
log.file <- file(description = paste0(PATH, "GA", GA.id, "_", CROP, ".log"),
                 open        = "wb",
                 encoding    = "UTF-8")
write_log(paste0("INITIALIZING GA", GA.id))

MU       <- 70L
MAX.SIZE <- Inf
MAX.GEN  <- 100L
P.RECOMB <- 0.7
P.MUT    <- 0.2
MAX.DEPTH <- 0.9 #1.1 #2.5 #
NUM.DEPTHS <- ((MAX.DEPTH - 0.1) / 0.2) + 1
ROUND <- FALSE

##### INPUTS #####
write_log("\nREADING INPUT FILES")
PARAMS <- readr::read_csv(paste0(input.path, "crop_water_calibration_parameters_", CROP, ".csv"), col_types = readr::cols()) %>%
  dplyr::filter(calibrate == TRUE) %>%
  dplyr::select(-calibrate)
INITIAL.POP <- ecr::genReal(MU,
                            n.dim = nrow(PARAMS),
                            lower = PARAMS$param.min,
                            upper = PARAMS$param.max)

##### GA CONTROLS #####
write_log("\nSETTING GA CONTROLS")
control <- ecr::initECRControl(fitness.fun  = crop_water_GA_fitness,
                               n.objectives = as.integer(NUM.DEPTHS + 1),
                               minimize     = TRUE) %>%
  ecr::registerECROperator(slot = "mutate",
                           fun  = ecr::setup(ecr::mutPolynomial,
                                             eta   = 25,
                                             p     = P.MUT,
                                             lower = PARAMS$param.min,
                                             upper = PARAMS$param.max)) %>%
  ecr::registerECROperator(slot = "recombine",
                           fun  = ecr::setup(ecr::recSBX,
                                             eta   = 15,
                                             p     = P.RECOMB,
                                             lower = PARAMS$param.min,
                                             upper = PARAMS$param.max)) %>%
  ecr::registerECROperator(slot = "selectForMating",
                           fun  = ecr::selSimple) # %>% #selGreedy
# ecr::registerECROperator(slot = "selectForSurvival",
#                          fun  = ecr::selNondom)

##### RUN STICS REFERENCE SIMULATION #####
write_log("\nRUNNING STICS REFERNECE SIMULATION")
dum <- file.copy(from      = paste0(template.path, "stics"),
                 to        = output.path,
                 recursive = TRUE,
                 overwrite = TRUE)

write_line <- function(x) cat(x, file = batch.script, sep = "\n", append = TRUE)

batch.script <- file(description = paste0(output.path, "stics.sh"),
                     open        = "wb",
                     encoding    = "UTF-8")
cat("", file = batch.script, sep = "", append = FALSE)
write_line("#!/bin/sh")
write_line("#SBATCH --account=hisafe")
write_line("#SBATCH --partition=defq")
#write_line("#SBATCH --account=f_meso")
#write_line("#SBATCH --partition=fmuse1")
write_line("#SBATCH --mail-type=NONE")
write_line("#SBATCH --mail-user=wolzkevin@gmail.com")
write_line("module purge")
write_line("module load jre/jre.8_x64")
write_line("cd /nfs/work/hisafe/Capsis4")
write_line(paste0("sh capsis.sh -p script safe.pgms.ScriptGen ", output.path, "stics/stics.sim water_calib_defaults"))
close(batch.script)

system(paste0("sbatch ", output.path, "stics.sh"))

incomplete <- TRUE
while(incomplete) {
  Sys.sleep(1)
  complete <- check_completion(1, paste0(output.path, "stics"))
  if(complete) incomplete <- FALSE
}

stics.cells  <- read_profiles(paste0(output.path, "stics"), "cells",  hisafe = FALSE)
stics.voxels <- read_profiles(paste0(output.path, "stics"), "voxels", hisafe = FALSE)

lai.output <- stics.cells %>%
  dplyr::select(Day, Month, Year, JulianDay, lai)
write_delim(lai.output, paste0(output.path, "lai.obs"), delim = "\t", col_names = FALSE)

GROWTH.DATES <- stics.cells %>%
  dplyr::filter(lai > 0) %>%
  .$Date %>%
  unique()

STICS <- stics.voxels %>%
  dplyr::filter(Date %in% GROWTH.DATES) %>%
  dplyr::select(Date, z, cropWaterUptake) %>%
  dplyr::rename(stics = cropWaterUptake)

##### RUN OLD WINNER SIMULATION #####
write_log("\nRUNNING OLD WINNER SIMULATION")
dum <- file.copy(from      = paste0(template.path, "hisafe-old"),
                 to        = output.path,
                 recursive = TRUE,
                 overwrite = TRUE)
dum <- file.copy(from = paste0(output.path, "lai.obs"),
                 to   = paste0(output.path, "hisafe-old"))

write_line <- function(x) cat(x, file = batch.script, sep = "\n", append = TRUE)

batch.script <- file(description = paste0(output.path, "hisafe-old.sh"),
                     open        = "wb",
                     encoding    = "UTF-8")
cat("", file = batch.script, sep = "", append = FALSE)
write_line("#!/bin/sh")
write_line("#SBATCH --account=hisafe")
write_line("#SBATCH --partition=defq")
#write_line("#SBATCH --account=f_meso")
#write_line("#SBATCH --partition=fmuse1")
write_line("#SBATCH --mail-type=NONE")
write_line("#SBATCH --mail-user=wolzkevin@gmail.com")
write_line("module purge")
write_line("module load jre/jre.8_x64")
write_line("cd /nfs/work/hisafe/Capsis4")
write_line(paste0("sh capsis.sh -p script safe.pgms.ScriptGen ", output.path, "hisafe-old/hisafe-old.sim water_calib_defaults"))
close(batch.script)

system(paste0("sbatch ", output.path, "hisafe-old.sh"))

incomplete <- TRUE
while(incomplete) {
  Sys.sleep(1)
  complete <- check_completion(1, paste0(output.path, "hisafe-old"))
  if(complete) incomplete <- FALSE
}

old.winner.voxels <- read_profiles(paste0(output.path, "hisafe-old"), "voxelsOptim", hisafe = FALSE)

OLD.WINNER <- old.winner.voxels %>%
  dplyr::filter(Date %in% GROWTH.DATES) %>%
  dplyr::select(Date, z, cropWaterUptake) %>%
  dplyr::rename(old.winner = cropWaterUptake)

##### EVALUATE INITIAL POPULATION #####
write_log("\nEVALUATING INITIAL POPULATION")
GEN <- 0
gen.folder <- paste0(gen.path, "gen", GEN, "/")
population <- INITIAL.POP
if(ROUND) population <- purrr::map(population, round_params)
fitness <- crop_water_GA_fitness(population)

initial.fitness.summary <- summarize_fitness(fitness)
readr::write_csv(initial.fitness.summary, paste0(output.path, "pop_fitness_summary.csv"))

##### INITIALIZE PARETO ARCHIVE #####
write_log("\nINITIALIZING PARETO ARCHIVE")
truncateByHVContr <- function(inds, fitness, max.size, ...) {
  hvcs     <- ecr::computeHVContr(fitness, ...)
  hvcs.ord <- order(hvcs, decreasing = TRUE)
  return(list(individuals = inds[hvcs.ord[seq_len(max.size)]],
              fitness     = fitness[, hvcs.ord[seq_len(max.size)], drop = FALSE]))
}
archive <- ecr::initParetoArchive(control   = control,
                                  max.size  = MAX.SIZE,
                                  trunc.fun = truncateByHVContr)

##### RUN GA #####
for(GEN in seq_len(MAX.GEN)) {
  write_log(paste("\nSTARTING GENERATION:", GEN))
  startGen <- proc.time()[3]

  gen.folder <- paste0(gen.path, "gen", GEN, "/")

  ## GET POP DATABASE
  write_log("- get population database")
  previous.pop.files <- paste0(output.path, "GA_data_", 0:(GEN-1), ".csv")
  pop.database <- purrr::map_df(previous.pop.files, readr::read_csv, col_types = cols())

  ## UPDATE PARETO ARCHIVE
  write_log("- update pareto archive")
  ecr::updateParetoArchive(archive = archive,
                           inds    = population,
                           fitness = fitness)
  fitness.archive <- ecr::getFront(archive)
  inds.archive    <- ecr::getIndividuals(archive)

  ## SUMMARIZE FITNESS OF PARETO SOLUTIONS
  write_log("- summarize pareto solutions")
  pareto.fitness.summary <- summarize_fitness(fitness.archive)
  readr::write_csv(pareto.fitness.summary, paste0(output.path, "pareto_fitness_summary.csv"), append = (GEN != 1))

  ## PLOT & WRITE OUT PARETO FRONT
  # data.files <- list.files(output.path, pattern = "GA_data_[0-9].", full.names = TRUE)
  # all.sim.data <- purrr::map_df(data.files, read_csv, col_types= cols())

  # which(dplyr::distinct(all.sim.data))
  # dplyr::distinct_at(PARAMS$param.name, .keep_all = TRUE)

  pareto.fitness <- fitness.archive %>%
    t() %>%
    as.data.frame()

  pareto.params <- inds.archive %>%
    unlist() %>%
    matrix(ncol = nrow(PARAMS), byrow = TRUE) %>%
    dplyr::as_tibble()
  names(pareto.params) <- PARAMS$param.name

  # pareto.params <- pareto.params %>%
  #   dplyr::left_join(all.sim.data, by = PARAMS$param.name)

  pareto.front <- pareto.params %>%
    dplyr::mutate(gen = GEN) %>%
    dplyr::bind_cols(pareto.fitness) %>%
    dplyr::select(gen, dplyr::everything())


  readr::write_csv(pareto.front, paste0(output.path, "pareto_front.csv"),  append = (GEN != 1))

  ## PLOT GA DIAGNOSTICS
  if(GEN > 1) {
    write_log("- plot GA diagnostics")
    fit.summary <- readr::read_csv(paste0(output.path, "pareto_fitness_summary.csv"), col_types = readr::cols())
    fit.trjectory <- ggplot2::ggplot(fit.summary, aes(x = gen)) +
      ggplot2::labs(x = "Generation", y = "Fitness") +
      ggplot2::facet_wrap(~metric, ncol = 1, scales = "free_y") +
      ggplot2::geom_ribbon(aes(ymin = min, ymax = max), fill = "grey90") +
      ggplot2::geom_line(aes(y = mean),   color = "black", linetype = "solid",  size = 1) +
      ggplot2::geom_line(aes(y = median), color = "black", linetype = "dashed", size = 1) +
      ggplot2::theme_bw()
    dum <- ggsave_fitmax(paste0(output.path, "GA_fitness_trajectory.pdf"), fit.trjectory)

    pareto.size <- fit.summary %>%
      dplyr::group_by(gen) %>%
      dplyr::summarize(n = mean(n))

    pareto.size.plot <- ggplot2::ggplot(pareto.size, aes(x = gen)) +
      ggplot2::labs(x = "Generation", y = "Size of Pareto archive") +
      ggplot2::geom_line(aes(y = n),   color = "black", linetype = "solid",  size = 1) +
      ggplot2::theme_bw()
    dum <- ggsave_fitmax(paste0(output.path, "GA_pareto_archive_size.pdf"), pareto.size.plot)

    pareto.pairs.plot <- GGally::ggpairs(pareto.params)
    dum <- ggplot2::ggsave(paste0(output.path, "GA_pareto_param_pairs_plot_", GEN, ".png"), pareto.pairs.plot)

    hist.data <- pareto.params %>%
      tidyr::gather(key = "param")

    sol.hist <- ggplot2::ggplot(hist.data, aes(x = value)) +
      ggplot2::geom_histogram(bins = 10) +
      ggplot2::facet_wrap(~param, scales = "free_x") +
      ggplot2::theme_bw()
    dum <- ggsave_fitmax(paste0(output.path, "GA_pareto_histogram_gen", GEN, ".pdf"), sol.hist)
  }

  if(GEN < MAX.GEN) {
    ## GENERATE NEW POPULATION
    write_log("- generate new population")
    population <- ecr::generateOffspring(control  = control,
                                         inds     = inds.archive,
                                         fitness  = fitness.archive,
                                         lambda   = MU,
                                         p.recomb = P.RECOMB,
                                         p.mut    = P.MUT)
    if(ROUND) population <- purrr::map(population, round_params)

    ## EVALUATE FITNESS OF NEW POPULATION
    write_log("- evaluate fitness of new population")
    fitness <- crop_water_GA_fitness(population)

    pop.fitness <- fitness %>%
      t() %>%
      as.data.frame()

    ## SUMMARIZE FITNESS OF POPULATION
    write_log("- summarize fitness of new population")
    pop.fitness.summary <- summarize_fitness(fitness)
    readr::write_csv(pop.fitness.summary, paste0(output.path, "pop_fitness_summary.csv"), append = TRUE)
  }

  ## SAVE PARETO ARCHIVE
  write_log("- save pareto archive")
  save(archive, file = paste0(output.path, "Pareto_Archive_gen", GEN, ".RData"))

  elapsedGen = round((proc.time()[3] - startGen))
  write_log(paste0("DONE WITH  GENERATION: ", GEN, " (", elapsedGen, " seconds)"))
}

##### RUN BEST #####
write_log("\nRERUN BEST SOLUTIONS")
GEN <- "MINDIF"
gen.folder <- paste0(gen.path, "gen", GEN, "/")
min.dif.id <- which.min(fitness.archive[nrow(fitness.archive),])
min.dif.params <- list(inds.archive[[min.dif.id]])
min.dif.fitness <- crop_water_GA_fitness(min.dif.params)
min.dif.fitness.summary <- summarize_fitness(min.dif.fitness)
readr::write_csv(min.dif.fitness.summary, paste0(output.path, "pop_fitness_summary.csv"), append = TRUE)

GEN <- "MINVOL"
gen.folder <- paste0(gen.path, "gen", GEN, "/")
min.vol.id <- which.min(apply(fitness.archive, 2, prod))
min.vol.params <- list(inds.archive[[min.vol.id]])
min.vol.fitness <- crop_water_GA_fitness(min.vol.params)
min.vol.fitness.summary <- summarize_fitness(min.vol.fitness)
readr::write_csv(min.vol.fitness.summary, paste0(output.path, "pop_fitness_summary.csv"), append = TRUE)

write_log(paste0("DONE WITH GA", GA.id))