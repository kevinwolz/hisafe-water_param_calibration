### WATER CALIBRATION ALGORITHM - PREP
### Author: Kevin J. Wolz

library(hisafer)

GA.id <- 1
YEARS <- 2007:2008
REFERENCE.PROFILES <- c("cells", "voxels", "plot", "cellsDetail")

for(CROP in c("durum-wheat-restinclieres", "weed-restinclieres", "rape", "winter-pea")) {

  PARAMS <- readr::read_csv(paste0(input.path, "crop_water_calibration_parameters_", CROP, ".csv"), col_types = readr::cols()) %>%
    dplyr::filter(calibrate == TRUE) %>%
    dplyr::select(-calibrate)


  PREP.PATH <- paste0("./water_calib_GA", GA.id, "/templates/", CROP)

  common.params <- list(nbSimulations        = length(YEARS),
                        waterTable           = 0,
                        simulationYearStart  = YEARS[1],
                        mainCropSpecies      = paste0(CROP, ".plt"),
                        mainCropItk          = paste0(CROP, ".tec"),
                        interCropSpecies     = paste0(CROP, ".plt"),
                        interCropItk         = paste0(CROP, ".tec"),
                        layers               = layer_params(template = "monocrop",
                                                            thick = c(0.4, 0.4, 0.6, 0.6)))

  ##### STICS REFERENCE SIMULATION #####
  ref.hip <- define_hisafe(path                 = PREP.PATH,
                           template             = "monocrop",
                           profiles             = REFERENCE.PROFILES,
                           SimulationName       = "stics",
                           sticsWaterExtraction = 1,
                           bulk.pass            = common.params)
  build_hisafe(ref.hip, plot.scene = FALSE, summary.files = FALSE, files = c("sim", "plt", "par"))

  ##### TEMPLATE HISAFE SIMULATION #####
  hip <- define_hisafe(path                 = PREP.PATH,
                       template             = "monocrop",
                       profiles             = "voxelsOptim",
                       SimulationName       = "hisafe",
                       sticsWaterExtraction = 0,
                       laiFileName          = "lai.obs",
                       bulk.pass            = common.params)
  build_hisafe(hip, plot.scene = FALSE, summary.files = FALSE, files = c("sim", "plt", "par"))

  ##### OLD HISAFE SIMULATION #####
  old.params <- as.list(PARAMS$old.winner)
  names(old.params) <- PARAMS$param.name

  hip <- define_hisafe(path                 = PREP.PATH,
                       template             = "monocrop",
                       profiles             = "voxelsOptim",
                       SimulationName       = "hisafe-old",
                       sticsWaterExtraction = 0,
                       laiFileName          = "lai.obs",
                       bulk.pass            = c(common.params, old.params))
  build_hisafe(hip, plot.scene = FALSE, summary.files = FALSE, files = c("sim", "plt", "par"))
}