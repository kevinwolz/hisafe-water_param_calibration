### WATER CALIBRATION ALGORITHM - PREP
### Author: Kevin J. Wolz

library(hisafer)
library(tidyverse)

GA.id <- 1
YEARS <- 2003:2005
REFERENCE.PROFILES <- c("cells", "voxels", "plot", "cellsDetail", "climate")
CROPS <- c("durum-wheat-restinclieres", "weed-restinclieres", "rape", "winter-pea")
input.path <- "./input/"

for(CROP in CROPS) {

  PREP.PATH <- paste0("./output", GA.id, "/templates/", CROP)
  dir.create(PREP.PATH, showWarnings = FALSE, recursive = TRUE)

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
  PARAMS <- readr::read_csv(paste0(input.path, "crop_water_calibration_parameters_", CROP, ".csv"), col_types = readr::cols()) %>%
    dplyr::select(-calibrate)

  old.params <- PARAMS %>%
    .$old.winner %>%
    as.list()

  names(old.params) <- PARAMS$param.name

  hip <- define_hisafe(path                 = PREP.PATH,
                       template             = "monocrop",
                       profiles             = REFERENCE.PROFILES,
                       SimulationName       = "hisafe-old",
                       sticsWaterExtraction = 0,
                       laiFileName          = "lai.obs",
                       bulk.pass            = c(common.params, old.params))
  build_hisafe(hip, plot.scene = FALSE, summary.files = FALSE, files = c("sim", "plt", "par"))
}

##### DEFAULTS #####
hip <- define_hisafe(path                 = ".",
                     template             = "monocrop",
                     profiles             = c(REFERENCE.PROFILES, "voxelsOptim"),
                     SimulationName       = "water_calib_defaults",
                     bulk.pass            =  list(waterTable           = 0,
                                                  nbSimulations        = length(CROPS),
                                                  mainCropSpecies      = list(paste0(CROPS, ".plt")), # just used to get the .plt files into the defaults
                                                  mainCropItk          = list(paste0(CROPS, ".tec")), # just used to get the .tec files into the defaults
                                                  layers               = layer_params(template = "monocrop",
                                                                                      thick = c(0.4, 0.4, 0.6, 0.6))))
build_hisafe(hip, plot.scene = FALSE, summary.files = FALSE, files = c("pld", "wth", "pro", "tec"))
