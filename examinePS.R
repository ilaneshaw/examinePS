## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "examinePS",
  description = "",
  keywords = "",
  authors = structure(list(list(given = c(""), family = "", role = c("aut", "cre"), email = "", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(examinePS = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("NEWS.md", "README.md", "examinePS.Rmd"),
  reqdPkgs = list(
    "PredictiveEcology/SpaDES.core@development (>= 2.0.2.9000)", "ggplot2", "sf", "data.table", "terra",
    "LandR", "googledrive", "plotrix", "ggpubr", "diptest", "nortest", "dplyr", "tidyverse", "reshape2", "gt", "gtExtras"
  ),
  parameters = bindrows(
    # defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(
      ".plots", "character", "screen", NA, NA,
      "Used by Plots function, which can be optionally used here"
    ),
    defineParameter(
      ".plotInitialTime", "numeric", start(sim), NA, NA,
      "Describes the simulation time at which the first plot event should occur."
    ),
    defineParameter(
      "doPredsInitialTime", "numeric", start(sim), NA, NA,
      "Describes the simulation time at which the first plot event should occur."
    ),
    defineParameter(
      "doPredsInterval", "numeric", NA, NA, NA,
      "Describes the simulation time interval between getPreds events."
    ),
    defineParameter(
      ".plotInterval", "numeric", NA, NA, NA,
      "Describes the simulation time interval between plot events."
    ),
    defineParameter(
      ".saveInitialTime", "numeric", NA, NA, NA,
      "Describes the simulation time at which the first save event should occur."
    ),
    defineParameter(
      ".saveInterval", "numeric", 1, NA, NA,
      "This describes the simulation time interval between save events."
    ),
    defineParameter(
      ".studyAreaName", "character", NA, NA, NA,
      "Human-readable name for the study area used - e.g., a hash of the study",
      "area obtained using `reproducible::studyAreaName()`"
    ),
    defineParameter(
      "only1DPS", "logical", FALSE, NA, NA,
      "do smoothing by cover class only (1D)? if FALSE smoothing will be done by forest type and age class where possible"
    ),
    defineParameter(
      "maxAgeClass", "numeric", 17, NA, NA,
      "what the oldest age class will be (everything older will be included in this class)"
    ),
    defineParameter(
      "ageGrouping", "numeric", 10, NA, NA,
      "how many years included per age class"
    ),
    defineParameter(
      "spList", "character", NA, NA, NA,
      "a list of species in the format of 4-letter codes"
    ),
    defineParameter(
      "rasterToMatchLocation", "character", NA, NA, NA,
      "the file location of the rasterToMatch"
    ),
    defineParameter(
      "rasterToMatchName", "character", NA, NA, NA,
      "the name of the rasterToMatch file"
    ),
    defineParameter(
      "studyAreaLocation", "character", NA, NA, NA,
      "the file location of the studyArea"
    ),
    defineParameter(
      "nameBCR", "character", NA, NA, NA,
      "the BAM regional model BCR region that the studyArea is located in"
    ),
    defineParameter(
      "nameForClassRas", "character", NA, NA, NA,
      "the file name of the forest class raster"
    ),
    defineParameter(
      "locationForClass", "character", NA, NA, NA,
      "the location of the forest class raster"
    ),
    defineParameter(
      "nameLandClassRas", "character", NA, NA, NA,
      "the file name of the non forest raster"
    ),
    defineParameter(
      "locationLandClass", "character", NA, NA, NA,
      "the location of the non forest raster"
    ),
    defineParameter(
      "nameAgeRas", "character", NA, NA, NA,
      "the file name of the age raster"
    ),
    defineParameter(
      "locationAge", "character", NA, NA, NA,
      "the location of the age raster"
    ),
    defineParameter(
      "locationSpRaster", "character", NA, NA, NA,
      "the location of the sp density rasters"
    ),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(
      ".seed", "list", list(), NA, NA,
      "Named list of seeds to use for each event (names)."
    ),
    defineParameter(
      ".useCache", "logical", FALSE, NA, NA,
      "Should caching of events or module be used?"
    ),
    defineParameter(
      "min2DStatsSample", "numeric", 100, 2, NA,
      "exclude any classes from 2D stats table that have a sample size smaller than minStatsSample"
    )
  ),
  inputObjects = bindrows(
    # expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    # expectsInput("rasterToMatch", "SpatRaster", desc = "A raster used to determine projection of other spatial objects. Must cover all of the region covered by the studyArea"),
    # expectsInput("studyArea", "SpatVector", desc = "Polygon to use as the study area."),
    # expectsInput(objectName = "forClassRaster", objectClass = "SpatRaster", desc = NA),
    # expectsInput(objectName = "landClassRaster", objectClass = "SpatRaster", desc = NA),
    # expectsInput(objectName = "ageRaster", objectClass = "SpatRaster", desc = NA),
    # expectsInput(objectName = "spRasters", objectClass = NA, desc = NA),
    # expectsInput(objectName = "spDatasets", objectClass = NA, desc = NA),
    # expectsInput(objectName = "spPreds", objectClass = NA, desc = NA),
    # expectsInput(objectName = "statsGBM", objectClass = NA, desc = NA),
    # expectsInput(objectName = "for1DMaps", objectClass = NA, desc = NA),
    # expectsInput(objectName = "for1DAndLc1DMaps", objectClass = NA, desc = NA),
    # expectsInput(objectName = "for2DAndLc1DMaps", objectClass = NA, desc = NA),
    # expectsInput(objectName = "for2DMaps", objectClass = NA, desc = NA)
  ),
  outputObjects = bindrows(
    # createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = NA, objectClass = NA, desc = NA),
    createsOutput(objectName = "assumpTab1D", objectClass = NA, desc = NA),
    createsOutput(objectName = "assumptionsByClass1D", objectClass = NA, desc = NA),
    createsOutput(objectName = "assumptionsBySp1D", objectClass = NA, desc = NA),
    createsOutput(objectName = "for1DRes", objectClass = NA, desc = NA),
    createsOutput(objectName = "for1DAndLc1DRes", objectClass = NA, desc = NA),
    createsOutput(objectName = "spStats2D", objectClass = NA, desc = NA),
    createsOutput(objectName = "assumptionsByClass2D", objectClass = NA, desc = NA),
    createsOutput(objectName = "assumptionsBySp2D", objectClass = NA, desc = NA),
    createsOutput(objectName = "spearmanStats", objectClass = NA, desc = NA),
    createsOutput(objectName = "residualStats", objectClass = NA, desc = NA),
    createsOutput(objectName = "for2DAndLc1DRes", objectClass = NA, desc = NA),
    createsOutput(objectName = "for2DRes", objectClass = NA, desc = NA)
  )
))

## event types
#   - type `init` is required for initialization

doEvent.examinePS <- function(sim, eventTime, eventType) {
  switch(eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim,
        eventTime = P(sim)$doPredsInitialTime,
        moduleName = "examinePS", eventType = "examine1D"
      )
      if (P(sim)$only1DPS == FALSE) {
        sim <- scheduleEvent(sim,
          eventTime = P(sim)$doPredsInitialTime,
          moduleName = "examinePS", eventType = "examine2D"
        )
        sim <- scheduleEvent(sim,
          eventTime = P(sim)$doPredsInitialTime,
          moduleName = "examinePS", eventType = "compare1D2D"
        )
      }
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "examinePS", "plot")
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "examinePS", "save")
    },
    plot = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      plotFun(sim) # example of a plotting function
      # schedule future event(s)

      # e.g.,
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "examinePS", "plot")

      # ! ----- STOP EDITING ----- ! #
    },
    save = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event

      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + P(sim)$.saveInterval, "examinePS", "save")

      # ! ----- STOP EDITING ----- ! #
    },
    examine1D = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- examine1D(sim)
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "examinePS", "templateEvent")
      sim <- scheduleEvent(sim,
        eventTime = time(sim) + P(sim)$doPredsInterval,
        moduleName = "examinePS", eventType = "examine1D"
      )
      # ! ----- STOP EDITING ----- ! #
    },
    examine2D = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- examine2D(sim)
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "examinePS", "templateEvent")
      sim <- scheduleEvent(sim,
        eventTime = time(sim) + P(sim)$doPredsInterval,
        moduleName = "examinePS", eventType = "examine2D"
      )
      # ! ----- STOP EDITING ----- ! #
    },
    compare1D2D = {
      # ! ----- EDIT BELOW ----- ! #
      # do stuff for this event
      sim <- compare1D2D(sim)
      # e.g., call your custom functions/methods here
      # you can define your own methods below this `doEvent` function

      # schedule future event(s)

      # e.g.,
      # sim <- scheduleEvent(sim, time(sim) + increment, "examinePS", "templateEvent")
      sim <- scheduleEvent(sim,
        eventTime = time(sim) + P(sim)$doPredsInterval,
        moduleName = "examinePS", eventType = "compare1D2D"
      )

      # ! ----- STOP EDITING ----- ! #
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
      "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'",
      sep = ""
    ))
  )
  return(invisible(sim))
}

## event functions
#   - keep event functions short and clean, modularize by calling subroutines from section below.

### template initialization
Init <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

### template for save events
Save <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  sim <- saveFiles(sim)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for plot events
plotFun <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # do stuff for this event
  # sampleData <- data.frame("TheSample" = sample(1:10, replace = TRUE))
  # Plots(sampleData, fn = ggplotFn) # needs ggplot2
  # browser()
  lapply(P(sim)$spList, FUN = function(sp) {
    print(sp)
    nmRas <- eval(parse(text = paste("sim$spRasters$", sp, sep = "")))
    mapRas1D <- eval(parse(text = paste("sim$for1DAndLc1DMaps$", sp, sep = "")))
    resRas1D <- eval(parse(text = paste("sim$for1DAndLc1DRes$", sp, sep = "")))
    if (P(sim)$only1DPS == FALSE) {
      mapRas2D <- eval(parse(text = paste("sim$for2DAndLc1DMaps$", sp, sep = "")))
      resRas2D <- eval(parse(text = paste("sim$for2DAndLc1DRes$", sp, sep = "")))
    }
    clearPlot()
    Plot(nmRas, title = paste("BAM Sp Density Prediction Model for ", sp, sep = ""), na.color = "white")
    Plot(mapRas1D, title = paste("1D Mapped Predictions for ", sp, sep = ""), na.color = "white")
    Plot(resRas1D, title = paste("1D Residuals for ", sp, sep = ""), na.color = "white")

    if (P(sim)$only1DPS == FALSE) {
      clearPlot()
      Plot(nmRas, title = paste("BAM Sp Density Prediction Model for ", sp, sep = ""), na.color = "white")
      Plot(mapRas2D, title = paste("2D Mapped Predictions for ", sp, sep = ""), na.color = "white")
      Plot(resRas2D, title = paste("2D Residuals for ", sp, sep = ""), na.color = "white")
    }
  })


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### template for your event1
examine1D <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test
  # get summary of assumptions/stats for 1D

  print("get assumptions summaries")
  # make single dataframe
  spPreds1DSingleFrame <- rbindlist(sim$spPreds1D)

  # Get table of binary normality and unimodality
  # If p value is less than or equal to 0.05 it fails the test and is not considered normal/unimodal
  # fail = 0, pass = 1
  assumpTab1D <- spPreds1DSingleFrame[, c(1, 2, 9, 11, 12)]
  assumpTab1D$normal <- NA
  assumpTab1D$normal[assumpTab1D$normality_p < 0.05] <- 0
  assumpTab1D$normal[assumpTab1D$normality_p == 0.05] <- 0
  assumpTab1D$normal[assumpTab1D$normality_p > 0.05] <- 1


  assumpTab1D$unimodal <- NA
  assumpTab1D$unimodal[assumpTab1D$unimodality_p < 0.05] <- 0
  assumpTab1D$unimodal[assumpTab1D$unimodality_p == 0.05] <- 0
  assumpTab1D$unimodal[assumpTab1D$unimodality_p > 0.05] <- 1
  assumpTab1D <- assumpTab1D[, c(1, 2, 5, 6, 7)]


  # we make assumption that if there is an NA, it is not normal/unimodal
  assumpTab1D$normal[is.na(assumpTab1D$normal) == TRUE] <- 0
  assumpTab1D$unimodal[is.na(assumpTab1D$unimodal) == TRUE] <- 0


  # save tab for furture graphs comparing with 2D assumptions
  write.csv(sim$assumpTab1D, file = file.path(outputFolderSpPreds, "assumpTab1D.csv"))

  # get table of prop of sp with p values under 0.05 per class
  print("get assumptions by class 1D")
  sim$assumptionsByClass1D <- assumpTab1D[order(landForClass)][, list(
    noSps = .N,
    propSpsNormal = mean(normal),
    propSpsUnimodal = mean(unimodal),
    smoothingType = "1D"
  ),
  by = landForClass
  ]
  write.csv(sim$assumptionsByClass1D, file = file.path(outputFolderSpPreds, "assumptionsByClass1D.csv"))

  # get table of sp giving prop of classes with p values under 0.05
  print("get assumptions by sp 1D")
  sim$assumptionsBySp1D <- assumpTab1D[order(species)][, list(
    noClasses = .N,
    propClassesNormal = mean(normal),
    propClassesUnimodal = mean(unimodal),
    smoothingType = "1D"
  ),
  by = species
  ]
  write.csv(sim$assumptionsBySp1D, file = file.path(outputFolderSpPreds, "assumptionsBySp1D.csv"))


  ### Make residual rasters of composite 1D predictions for forClassraster areas and 1D predictions for landClassRaster areas
  print("Make for1DAndLc1DRes ")
  sim$for1DAndLc1DRes <- lapply(X = P(sim)$spList, FUN = function(sp) {
    print(sp)
    NM <- eval(parse(text = paste("sim$spRasters$", sp, sep = "")))
    PS <- eval(parse(text = paste("sim$for1DAndLc1DMaps$", sp, sep = "")))
    res <- NM - PS

    names(res) <- paste(sp)
    #
    # clearPlot()
    # Plot(res, na.colour = "blue", title = paste(sp,  " 1DPS residuals", sep = ""))
    #

    print(paste(sp, " for 1D and lc 1D res raster complete"))
    return(res)
  })

  names(sim$for1DAndLc1DRes) <- P(sim)$spList


  # save residual rasters
  lapply(X = P(sim)$spList, FUN = function(sp) {
    raster <- eval(parse(text = paste("sim$for1DAndLc1DRes$", sp, sep = "")))
    names(raster) <- paste(sp)
    terra::writeRaster(
      x = raster,
      filename = file.path(outputFolderSpPredsRasters, paste(sp, "-for1DAndLc1DRes", sep = "")),
      filetype = "GTiff",
      gdal = "COMPRESS=NONE",
      overwrite = TRUE
    )
  })


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}


### template for your event2
examine2D <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event2Test1 <- " this is test for event 2. " # for dummy unit test
  # sim$event2Test2 <- 777  # for dummy unit test
  # create spStats2D (MODULE OUTPUT), a list of lists giving a spStats2D data table and vector classesNotPresent for each sp species on the spList.
  print("get spStats2D")

  sim$spStats2D <- lapply(X = P(sim)$spList, FUN = function(sp) {
    print(sp)
    forestedDT <- as.data.table(eval(parse(text = paste("sim$spDatasets$", sp, sep = ""))))
    # separate out data table rows that are forested from the raw spDataset
    forestedDT <- forestedDT[FoLRaster == "forClass"]
    forestedDT <- droplevels(forestedDT)

    # get rid of any rows with NA for age
    forestedDT <- na.omit(forestedDT, cols = "age")

    # get age classes for each row using the ageClassDefs table
    ageClass <- forestedDT[, age]
    ageClass <- as.data.table(ageClass)
    ageClass[] <- lapply(ageClass, function(x) sim$spPreds$ageClassDefs$ageClasses[match(x, sim$spPreds$ageClassDefs$allAges)])
    # add the ageClass to the forestedDT
    spDataNew <- cbind(forestedDT, ageClass)

    # create new column, landAgeClass, giving the uniqueLandClass and ageClass combined
    spDataNew <- as.data.table(unite(spDataNew, landAgeClass, c(landForClass, ageClass), sep = ".", remove = FALSE))

    # produce data table of statistics on the sp Data based on the 2D bins
    singleSpStats2D <- spDataNew[
      order(landAgeClass) # order the rows by the land cover class
    ][, list(
      classCount = .N, # get the number of cells each cover class
      meanSpDensity = mean(spDensity), # get mean sp density
      medianSpDensity = median(spDensity),
      varSpDensity = var(spDensity) * (.N - 1) / .N, # get the variance for sp density for each class
      seSpDensity = std.error(spDensity), # get the standard error for sp density for each  class
      normality.p = tryCatch(ad.test(spDensity)$p.value, error = function(cond) {
        return(NaN)
      }), # ifelse(mean(spDensity) > 0, tryCatch(ad.test(spDensity)$p.value,error = function(cond){return(NA)}), NA),
      unimodality.p = dip.test(spDensity)$p.value,
      species = sp
    ),
    by = list(landAgeClass)
    ]

    # exclude any classes from table that have a sample size smaller than minStatsSample, a parameter
    singleSpStats2D <- subset(singleSpStats2D, classCount > P(sim)$min2DStatsSample)

    # #### get list of missing classes
    # landAgeClassesPresent <- unique(singleSpStats2D$landAgeClass) #get classes that are represented in spStats2D
    #
    # #get all classes possible
    # landClasses <- rep(unique(spDataNew$uniqueClasses), times = maxAgeClass)
    # landClasses <- as.data.table(landClasses)
    # ageClassReps <- rep(1:maxAgeClass, times = length(unique(spDataNew$uniqueClasses)))
    # ageClassReps <- as.data.table(ageClassReps)
    # allPossibleClasses <- cbind(landClasses, ageClassReps)
    # allPossibleClasses <-  unite(allPossibleClasses, allPossibleClasses, c(landClasses, ageClassReps), sep= ".", remove=FALSE)
    #
    # #get the classes not present in spStats2D
    # classesNotPresent <- setdiff(allPossibleClasses$allPossibleClasses, landAgeClassesPresent)
    #
    # ##make list object of all stats outputs
    # spStatsList2D <- list(singleSpStats2D, classesNotPresent)
    # names(spStatsList2D) <- c("spStats2D", "classesNotPresent")

    # return(spStatsList2D)

    return(singleSpStats2D)
  })

  names(sim$spStats2D) <- P(sim)$spList

  # get 2D stats Summary
  print("get stats summary 2D")
  # make single dataframe of 2D stats
  spStats2DSingleFrame <- rbindlist(sim$spStats2D)

  # Get table of binary normality and unimodality
  assumpTab2D <- spStats2DSingleFrame[, c(1, 7, 8, 9)]
  assumpTab2D$normal <- NA
  assumpTab2D$normal[assumpTab2D$normality > 0.05] <- 0
  assumpTab2D$normal[assumpTab2D$normality == 0.05] <- 1
  assumpTab2D$normal[assumpTab2D$normality < 0.05] <- 1


  assumpTab2D$unimodal <- NA
  assumpTab2D$unimodal[assumpTab2D$unimodality > 0.05] <- 0
  assumpTab2D$unimodal[assumpTab2D$unimodality == 0.05] <- 1
  assumpTab2D$unimodal[assumpTab2D$unimodality < 0.05] <- 1
  assumpTab2D <- assumpTab2D[, c(1, 4, 5, 6)]
  assumpTab2D
  write.csv(assumpTab2D, file = file.path(outputFolderSpPreds, "assumpTab2D.csv"))

  # get table of prop of sp with p values under 0.05 per class
  sim$assumptionsByClass2D <- assumpTab2D[order(landAgeClass)][, list(
    noSps = .N,
    propSpsNormal = mean(normal),
    propSpsUnimodal = mean(unimodal),
    smoothingType = "2D"
  ),
  by = landAgeClass
  ]
  write.csv(sim$assumptionsByClass2D, file = file.path(outputFolderSpPreds, "assumptionsByClass2D.csv"))

  # get table of sp giving prop of classes with p values under 0.05
  sim$assumptionsBySp2D <- assumpTab2D[order(species)][, list(
    noClasses = .N,
    propClassesNormal = mean(normal),
    propClassesUnimodal = mean(unimodal),
    smoothingType = "2D"
  ),
  by = species
  ]
  write.csv(sim$assumptionsBySp2D, file = file.path(outputFolderSpPreds, "assumptionsBySp2D.csv"))


  # Get residual rasters

  # Make residual rasters of composite 2D predictions for forClassraster areas and 1D predictions for landClassRaster areas
  print("make for2DAndLc1DRes")
  sim$for2DAndLc1DRes <- lapply(X = P(sim)$spList, FUN = function(sp) {
    print(sp)

    NM <- eval(parse(text = paste("sim$spRasters$", sp, sep = "")))
    PS <- eval(parse(text = paste("sim$for2DAndLc1DMaps$", sp, sep = "")))
    res <- NM - PS

    names(res) <- paste(sp)
    # clearPlot()
    # Plot(res, na.colour = "blue", title = paste(sp,  " 2D smoothing residuals", sep = ""))
    #
    print(paste(sp, " for 2D and lc 1D res raster complete"))
    return(res)
  })

  names(sim$for2DAndLc1DRes) <- P(sim)$spList

  # save rasters
  lapply(X = P(sim)$spList, FUN = function(sp) {
    raster <- eval(parse(text = paste("sim$for2DAndLc1DRes$", sp, sep = "")))
    names(raster) <- paste(sp)
    terra::writeRaster(
      x = raster,
      filename = file.path(outputFolderSpPredsRasters, paste(sp, "-for2DAndLc1DRes", sep = "")),
      filetype = "GTiff",
      gdal = "COMPRESS=NONE",
      overwrite = TRUE
    )
  })

  ### ANALYSIS

  # calculate spearman stats
  print("get spearman stats")
  spearmanStats <- lapply(X = P(sim)$spList, FUN = function(sp) {
    print(sp)

    nmRas <- eval(parse(text = paste("sim$spRasters$", sp, sep = "")))
    map1D <- eval(parse(text = paste("sim$for1DAndLc1DMaps$", sp, sep = "")))
    map2D <- eval(parse(text = paste("sim$for2DAndLc1DMaps$", sp, sep = "")))


    valsNM <- as.data.table(terra::values(nmRas, dataframe = FALSE))
    vals1DMap <- as.data.table(terra::values(map1D, dataframe = FALSE))
    vals2DMap <- as.data.table(terra::values(map2D, dataframe = FALSE))
    valsMaps <- cbind(valsNM, vals1DMap, vals2DMap)
    valsMaps <- na.omit(valsMaps)
    colnames(valsMaps) <- c("valsNM", "vals1DMap", "vals2DMap")
    # valsMaps <- as.data.table(valsMaps)
    head(valsMaps)
    # Check normality assumption
    # Shapiro-Wilk normality test for all data
    # ad.test(valsMaps$valsNM) # => p = 0.1229
    # ad.test(valsMaps$vals1DMap) # => p = 0.09
    # ad.test(valsMaps$vals2DMap)

    # library("ggpubr")
    # ggqqplot(valsMaps$valsNM, ylab = "National Model Prediction")
    # ggqqplot(valsMaps$vals1DMap, ylab = "1D Map Prediction")
    # ggqqplot(valsMaps$vals2DMap, ylab = "2D Map Prediction")

    spearman1D <- cor(valsMaps$valsNM, valsMaps$vals1DMap, method = "spearman")
    print(spearman1D)
    spearman2D <- cor(valsMaps$valsNM, valsMaps$vals2DMap, method = "spearman")
    print(spearman2D)

    spearmanStats <- matrix(c(spearman1D, spearman2D), ncol = 2, byrow = TRUE)
    colnames(spearmanStats) <- c("spearman1D", "spearman2D")
    row.names(spearmanStats) <- sp

    return(spearmanStats)
  })

  sim$spearmanStats <- do.call(rbind, spearmanStats)

  fileName <- "spearmanStats.csv"
  write.csv(sim$spearmanStats, file = file.path(outputFolderSpPreds, fileName))

  head(sim$spearmanStats)
  browser()

  # get tables of residuals
  print("get resTabs")
  sim$resTabs <- lapply(X = spList, FUN = function(sp) {
    ras1D <- eval(parse(text = paste("sim$for1DAndLc1DRes$", sp, sep = "")))
    ras2D <- eval(parse(text = paste("sim$for2DAndLc1DRes$", sp, sep = "")))

    resVals1D <- as.data.table(terra::values(ras1D, dataframe = FALSE))
    resVals1D <- setnames(resVals1D, "resVals")
    resVals1D <- na.omit(resVals1D)
    res1DLab <- rep("res1D", nrow(resVals1D))
    resVals1D <- cbind(resVals1D, binningType = res1DLab)

    resVals2D <- as.data.table(terra::values(ras2D, dataframe = FALSE))
    resVals2D <- setnames(resVals2D, "resVals")
    resVals2D <- na.omit(resVals2D)
    res2DLab <- rep("res2D", nrow(resVals2D))
    resVals2D <- cbind(resVals2D, binningType = res2DLab)

    resVals <- rbind(resVals1D, resVals2D)
    species <- rep(paste(sp), nrow(resVals))
    resVals <- cbind(resVals, species = species)
    resVals[, absResVals := abs(resVals)]
    print(resVals)

    # save table
    # fileName <- paste(sp, "_fullResDataset.csv")
    # write.csv(resVals, file = file.path(outputFolderSpPreds, fileName))

    return(resVals)
  })
  names(sim$resTabs) <- P(sim)$spList

  # get residual stats
  print("get residual stats")
  residualStats <- lapply(X = P(sim)$spList, FUN = function(sp) {
    print(sp)

    nmRas <- eval(parse(text = paste("sim$spRasters$", sp, sep = "")))
    res1D <- eval(parse(text = paste("sim$for1DAndLc1DRes$", sp, sep = "")))
    res2D <- eval(parse(text = paste("sim$for2DAndLc1DRes$", sp, sep = "")))

    resTab <- eval(parse(text = paste("sim$resTabs$", sp, sep = "")))
    absResStats <- resTab[, list(
      medResAbs = median(absResVals),
      maxResAbs = max(absResVals)
    ),
    by = binningType
    ]

    m1D <- median(terra::values(res1D, dataframe = FALSE), na.rm = TRUE)
    m2D <- median(terra::values(res2D, dataframe = FALSE), na.rm = TRUE)
    sa1D <- terra::autocor(res1D, method = "moran")
    sa2D <- terra::autocor(res2D, method = "moran")
    m1DAbs <- absResStats[binningType == "res1D"]$medResAbs
    m2DAbs <- absResStats[binningType == "res2D"]$medResAbs
    max1DAbs <- absResStats[binningType == "res1D"]$maxResAbs
    max2DAbs <- absResStats[binningType == "res2D"]$maxResAbs
    mNM <- median(terra::values(nmRas, dataframe = FALSE), na.rm = TRUE)
    medAbsProp1D <- m1DAbs / mNM
    medAbsProp2D <- m2DAbs / mNM

    residualStats <- matrix(c(m1D, m2D, sa1D, sa2D, m1DAbs, m2DAbs, max1DAbs, max2DAbs, mNM, medAbsProp1D, medAbsProp2D), ncol = 11, byrow = TRUE)
    colnames(residualStats) <- c("median1DRes", "median2DRes", "autocor1DRes", "autocor2DRes", "median1DResAbs", "median2DResAbs", "max1DResAbs", "max2DResAbs", "nmMedian", "propOfNM1DMed", "propOfNM2DMed")
    row.names(residualStats) <- sp
    print(paste(sp, " calculation complete"))
    return(residualStats)
  })

  sim$residualStats <- do.call(rbind, residualStats)

  fileName <- "resStats.csv"
  write.csv(sim$residualStats, file = file.path(outputFolderSpPreds, fileName))

  head(sim$residualStats)


  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

compare1D2D <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  # THE NEXT TWO LINES ARE FOR DUMMY UNIT TESTS; CHANGE OR DELETE THEM.
  # sim$event1Test1 <- " this is test for event 1. " # for dummy unit test
  # sim$event1Test2 <- 999 # for dummy unit test

  browser()
  # COMPARE UNIMODALITY AND NORMALITY
  names(sim$assumptionsByClass1D)[names(sim$assumptionsByClass1D) == "landForClass"] <- "landAgeClass"
  sim$assumptionsByClass1D <- sim$assumptionsByClass1D[1:6]
  sim$assumptionsByClass1D <- droplevels(sim$assumptionsByClass1D)
  # assumptionsByClass1D_land <- assumptionsByClass1D[1:6]
  # assumptionsByClass1D_land  <- droplevels(assumptionsByClass1D_land)
  assumptionsByClass <- rbind(sim$assumptionsByClass1D, sim$assumptionsByClass2D)
  # assumptionsByClass <- assumptionsByClass[,2:5]
  assumptionsByClass$smoothingType <- as.factor(assumptionsByClass$smoothingType)
  assumptionsByClass$landAgeClass <- as.factor(assumptionsByClass$landAgeClass)

  # get summaries
  print("assumptions by class 1D summary")
  summary(sim$assumptionsByClass1D)
  # gt_plt_summary(sim$assumptionsByClass1D)
  sdUnimodality_1D <- sd(sim$assumptionsByClass1D$propSpsUnimodal)
  sdUnimodality_2D <- sd(sim$assumptionsByClass2D$propSpsUnimodal)
  sdNormality_1D <- sd(sim$assumptionsByClass1D$propSpsNormal)
  sdNormality_2D <- sd(sim$assumptionsByClass2D$propSpsNormal)
  sim$sdAssumptions <- data.table(sdUnimodality_1D, sdUnimodality_2D, sdNormality_1D, sdNormality_2D)
  print(sim$sdAssumptions)

  print("assumptions by class 2D summary")
  summary(sim$assumptionsByClass2D)
  # gt_plt_summary(sim$assumptionsByClass2D)

  # COMPARE 1D vs 2D  UNIMODALITY
  print("1D vs 2D unimodality")
  sim$unimodalityTTest <- t.test(propSpsUnimodal ~ smoothingType,
    data = assumptionsByClass
  )
  print(sim$unimodalityTTest)

  # # COMPARE 1D vs 2D  NORMALITY
  print("1D vs 2D normality")
  sim$normalityTTest <- tryCatch(t.test(propSpsNormal ~ smoothingType,
    data = assumptionsByClass
  ), error = function(cond) {
    return(NaN)
  })
  print(sim$normalityTTest)


  # COMPARE SIMILARITY BETWEEN NM AND PS RASTERS FOR 1D VS 2D (SPEARMAN RANK)
  sim$spearmanStats <- as.data.frame(sim$spearmanStats)
  print("spearman stats summary")
  print(summary(sim$spearmanStats))
  # gt_plt_summary(sim$spearmanStats)
  sdSpearman1D <- sd(sim$spearmanStats$spearman1D)
  sdSpearman2D <- sd(sim$spearmanStats$spearman2D)
  sim$sdSpearman <- data.table(sdSpearman1D, sdSpearman2D)
  print(sim$sdSpearman)

  # names(sim$spearmanStats) <- c("species", "spearman1D", "spearman2D")
  # spearmanStats$ageClasses <- as.numeric(spearmanStats$ageClasses)
  spearmanTab <- gather(data = sim$spearmanStats, key = "smoothingType", value = "spearmanStat")
  # spearmanStats <- unite(spearmanStats, dataset, c(binningType, ageClasses), remove=FALSE)
  spearmanTab <- na.omit(spearmanTab)
  spearmanTab <- as.data.table(spearmanTab)


  # 1D vs 2D correlation
  print("compare 1D vs 2D spearman")
  sim$corrTTest <- tryCatch(t.test(spearmanStat ~ smoothingType,
    data = spearmanTab
  ), error = function(cond) {
    return(NaN)
  })
  print(sim$corrTTest)


  # COMPARE SPATIAL AUTOCORRELATION OF RESIDUALS
  # names(sim$residualStats) <- c("species", "medianRes1D", "medianRes2D","moran1D", "moran2D")
  sim$residualStats <- data.table(sim$residualStats)
  moranStats <- sim$residualStats[, c(3:4)]
  print("moran stats summary")
  print(summary(moranStats))
  sdMoran1D <- sd(moranStats$autocor1DRes)
  sdMoran2D <- sd(moranStats$autocor2DRes)
  sim$sdMoran <- data.table(sdMoran1D, sdMoran2D)
  print(sim$sdMoran)
  # gt_plt_summary(moranStats)
  moranStats <- gather(data = moranStats, key = "smoothingType", value = "MoransI")
  moranStats <- na.omit(moranStats)
  moranStats <- as.data.table(moranStats)


  # 1D vs 2D spatial autocorrelation of residuals

  print("compare 1D vs 2D Moran")
  sim$saTTest <- tryCatch(t.test(MoransI ~ smoothingType,
    data = moranStats
  ), error = function(cond) {
    return(NaN)
  })
  print(sim$saTTest)

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  # cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #


  # get rasterToMatch
  if (!suppliedElsewhere("rasterToMatch", sim)) {
    print("get rasterTomatch from local drive")
    sim$rasterToMatch <- terra::rast(file.path(P(sim)$rasterToMatchLocation, P(sim)$rasterToMatchName))
  }

  # get studyArea shapefile
  if (!suppliedElsewhere("studyArea", sim)) {
    print("get studyArea shapefile from local drive")
    studyArea <- terra::vect(file.path(P(sim)$studyAreaLocation, P(sim)$.studyAreaName))

    # postProcess studyArea
    # Cache(projectRaster, raster, crs = crs(newRaster))
    message("cache studyArea")
    sim$studyArea <- reproducible::Cache(reproducible::postProcess, studyArea,
      # destinationPath = P(sim)$studyAreaLocation,
      # filename2 = "studyArea",
      useTerra = TRUE,
      fun = "terra::vect", # use the function vect
      targetCRS = crs(sim$rasterToMatch), # make crs same as rasterToMatch
      # overwrite = FALSE,
      verbose = TRUE
    )
  }

  # crop and mask rasterToMatch
  sim$rasterToMatch <- reproducible::Cache(terra::crop, sim$rasterToMatch, sim$studyArea)
  sim$rasterToMatch <- reproducible::Cache(terra::mask, sim$rasterToMatch, sim$studyArea)
  names(sim$rasterToMatch) <- "rasterToMatch"

  # get forest class raster
  if (!suppliedElsewhere("forClassRaster", sim)) {
    print("get forClassRaster from local Drive")
    sim$forClassRaster <- terra::rast(file.path(P(sim)$locationForClass, P(sim)$nameForClassRas))
    sim$forClassRaster <- reproducible::Cache(reproducible::postProcess, sim$forClassRaster,
      # destinationPath = downloadFolderForestClass,
      # use the function raster
      # targetCRS = crs(sim$rasterToMatch),
      fun = "terra::rast",
      useTerra = TRUE,
      # use the specified rasterToMatch to reproject to
      rasterToMatch = sim$rasterToMatch,
      studyArea = sim$studyArea,
      useCache = getOption("reproducible.useCache", TRUE),
      # overwrite = TRUE,
      verbose = TRUE
    )

    names(sim$forClassRaster) <- c("forClassRaster")

    sim$forClassRaster[sim$forClassRaster == 0] <- NA
  }

  # get land cover raster
  if (!suppliedElsewhere("landClassRaster", sim)) {
    print("get landClassRaster from Drive")
    sim$landClassRaster <- terra::rast(file.path(P(sim)$locationLandClass, P(sim)$nameLandClassRas))
    sim$landClassRaster <- reproducible::Cache(terra::crop, sim$landClassRaster, sim$studyArea)
    sim$landClassRaster <- reproducible::Cache(terra::mask, sim$landClassRaster, sim$studyArea)
    # sim$landClassRaster <- postProcess(sim$landClassRaster,
    #                                 #destinationPath = downloadFolderForestClass,
    #                                 #use the function raster
    #                                 #targetCRS = crs(sim$rasterToMatch),
    #                                 fun = "terra::rast",
    #                                 useTerra = TRUE,
    #                                 #use the specified rasterToMatch to reproject to
    #                                 rasterToMatch = sim$rasterToMatch,
    #                                 studyArea = sim$studyArea,
    #                                 useCache = getOption("reproducible.useCache", TRUE),
    #                                 #overwrite = TRUE,
    #                                 verbose = TRUE)
    #
    # landClassRaster[landClassRaster == 0] <- NA

    names(sim$landClassRaster) <- c("landClassRaster")
    print(sim$landClassRaster)
    Plot(sim$landClassRaster, na.color = "blue")
    print(terra::unique(sim$landClassRaster))
    sim$landClassRaster <- terra::mask(sim$landClassRaster, sim$forClassRaster, inverse = TRUE)

    print(sim$landClassRaster)
    Plot(sim$landClassRaster, na.color = "blue")
    print(terra::unique(sim$landClassRaster))
  }

  if (P(sim)$only1DPS == FALSE) {
    # get ageRaster
    if (!suppliedElsewhere("ageRaster", sim)) {
      print("get ageRaster from Drive")
      ageRaster <- terra::rast(file.path(P(sim)$locationAge, P(sim)$nameAgeRas))
      sim$ageRaster <- reproducible::Cache(reproducible::postProcess, ageRaster,
        # destinationPath = downloadFolderForestClass,
        # use the function raster
        useTerra = TRUE,
        fun = "terra::rast",
        # targetCRS = crs(sim$rasterToMatch),
        # use the specified rasterToMatch to reproject to
        rasterToMatch = sim$rasterToMatch,
        # studyArea = sim$studyArea,
        useCache = getOption("reproducible.useCache", TRUE),
        overwrite = FALSE,
        verbose = TRUE
      )

      names(sim$ageRaster) <- c("ageRaster")
    }
  }


  # get sp density rasters
  if (!suppliedElsewhere("spRasters", sim)) {
    P(sim)$spList <- sort(P(sim)$spList)
    ## for each item in turn from rastersForSplist the following function is applied:
    sim$spRasters <-
      lapply(
        X = P(sim)$spList,
        FUN = function(sp) {
          print(sp)
          # sp <- paste(sp, ".tif", sep = "")
          rasterFile <- paste(sp, "-meanBoot_60", sep = "")
          return(terra::rast(file.path(locationSpRas, rasterFile)))
        }
      )

    # get the species codes as names for the downloadedRasters object, rather than using the whole filepath
    # X <- lapply(sim$rastersForSpList, substr, 8, 11) #works for strings of the form "mosaic-XXXX-run3.tif"
    X <- lapply(sim$rastersForSpList, substr, 1, 4) # works for strings of the form "XXXX-meanBoot.tif"
    names(sim$spRasters) <- X

    sim$spRasters <- lapply(X = sim$spRasters, FUN = function(RasterLayer) {
      ## the function postProcesses the layer, cropping and masking it to a given study area and rasterToMatch, and saving it to a given destination path

      print(RasterLayer)
      proRaster <- reproducible::postProcess(RasterLayer,
        # studyArea = sim$studyArea,
        rasterToMatch = sim$rasterToMatch,
        useTerra = TRUE,
        fun = "terra::rast",
        # destinationPath = downloadFolderSp,
        # filename2 = paste(downloadFolderSp, "/", names(RasterLayer), ".tif", sep = ""),
        # overwrite = TRUE,
        verbose = TRUE
      )
      # clearPlot()
      # Plot(proRaster, na.color= "grey")
      return(proRaster)
    })

    names(sim$spRasters) <- P(sim)$spList
  }

  if (P(sim)$only1DPS == FALSE) {
    if (!suppliedElsewhere("ageClassRaster", sim)) {
      print("You need to make an ageClassRaster")
    }
  }

  if (!suppliedElsewhere("spDatasets", sim)) {
    print("You are missing your PS module outputs")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

ggplotFn <- function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes(TheSample)) +
    ggplot2::geom_histogram(...)
}

### add additional events as needed by copy/pasting from above
