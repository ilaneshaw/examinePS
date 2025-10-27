library("conflicted")
library("SpaDES.core")

options(spades.useRequire = FALSE)

## make a list of directory paths
inputsDir <- checkPath("../../inputs", create = TRUE)
outputsDir <- checkPath("../../outputs", create = TRUE)
downloadFolderArea <- checkPath(file.path(inputsDir, "studyArea/studyArea_AB_BCR6"), create = TRUE)
downloadFolderForestClass <- checkPath(file.path(inputsDir, "forestClassRasters"), create = TRUE)
downloadFolderBird <- checkPath(file.path(inputsDir, "birdRasterFiles"), create = TRUE)
outputFolderBirdPreds <- checkPath(file.path(outputsDir, "outputBirdPreds"), create = TRUE)
outputFolderBirdPredsRasters <- checkPath(file.path(outputsDir, "outputBirdPredsRasters"), create = TRUE)
setPaths(
  modulePath = file.path("../../modules"),
  cachePath = file.path("../../cache"),
  scratchPath = file.path("../../scratch"),
  inputPath = inputsDir,
  outputPath = outputsDir
)

simPaths <- getPaths()



# birdList <- sort(c("CAWA", "DEJU"))
# new AB study Area bird list (25% prob of occurrence in 1% of the area, only those with 32 replicates)
birdList <- sort(c(
  "ALFL", "BBWA", "BCCH", "BOCH", "BRCR", "CMWA", "COYE",
  "DEJU", "GCKI", "GRAJ", "LEFL", "MOWA", "OVEN", "PAWA",
  "RBNU", "RCKI", "REVI", "SWTH", "TEWA", "YRWA"
))

rasterToMatchLocation <- inputsDir
rasterToMatchName <- "LCC2005_V1_4a.tif"
studyAreaLocation <- downloadFolderArea
nameBCR <- "60"

.studyAreaName <- "studyArea_AB_BCR6.shp"

nameForClassRaster <- "vegTypesRas_AB_BCR6_2011"
folderUrlForClass <- downloadFolderForestClass

nameLandClassRaster <- "landCoverRas_AB_BCR6_2010"
folderUrlLandClass <- downloadFolderForestClass

nameAgeRaster <- "ageRas_AB_BCR6_2011"
folderUrlAge <- downloadFolderForestClass
folderUrlBirdRaster <- downloadFolderBird


simModules <- list("PS", "examinePS")

## Set simulation and module parameters
simTimes <- list(start = 1, end = 1, timeunit = "year")
simParams <- list(
  PS = list(
    doPredsInitialTime = 1,
    .plotInitialTime = 1,
    .saveInitialTime = 1,
    classOnly = FALSE,
    nTrees = 5000,
    ageGrouping = 10,
    maxAgeClass = 17,
    birdList = birdList,
    .studyAreaName = .studyAreaName,
    rasterToMatchLocation = rasterToMatchLocation,
    rasterToMatchName = rasterToMatchName,
    studyAreaLocation = studyAreaLocation,
    nameBCR = nameBCR,
    nameForClassRaster = nameForClassRaster,
    folderUrlForClass = folderUrlForClass,
    nameLandClassRaster = nameLandClassRaster,
    folderUrlLandClass = folderUrlLandClass,
    nameAgeRaster = nameAgeRaster,
    folderUrlAge = folderUrlAge
  ),
  examinePS = list(
    doPredsInitialTime = 1,
    .plotInitialTime = 1,
    .saveInitialTime = 1,
    fromDrive = FALSE,
    classOnly = FALSE,
    birdList = birdList
  )
)


## Simulation setup
mySim <- simInit(
  times = simTimes, params = simParams,
  modules = simModules, paths = simPaths
)

test <- spades(mySim)
