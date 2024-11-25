#' -----------------------------------------------------------------------------------------------------------------
#' Script: Helper Functions for Step 4 (S4)
#' Author: Caio Graco-Roza
#' Last Updated: 2024-11-24
#' -----------------------------------------------------------------------------------------------------------------
#' Description:
#' This script contains auxiliary functions for running Bayesian Bootstrap Generalized Dissimilarity Models (BBGDM) 
#' and handling associated data transformations. The functions assist in analyzing dissimilarity metrics and 
#' environmental predictors, including residual diagnostics, plotting, and effect size computations.
#' -----------------------------------------------------------------------------------------------------------------

require(Hmisc)
#' @title Pooled Variance Calculation
#' @description Calculates a pooled variance measure for two samples using a specified function.
#' @param x A numeric vector representing the first sample.
#' @param y A numeric vector representing the second sample.
#' @param FUN A function to calculate the variance measure (e.g., standard deviation).
#' @return A numeric value of the pooled variance measure.
#' @examples
#' pooled(c(1, 2, 3), c(4, 5, 6), sd)
pooled <- function(x, y, FUN) {
  nx <- length(x)
  ny <- length(y)
  sqrt(((nx - 1) * FUN(x) ^ 2 + (ny - 1) * FUN(y) ^ 2) / (nx + ny - 2))
}

#' @title High-Density Median
#' @description Calculates the high-density (HD) median of a numeric vector.
#' @param x A numeric vector.
#' @return A numeric value representing the HD median.
#' @examples
#' hdmedian(c(1, 2, 3, 4, 5))
hdmedian <- function(x) as.numeric(hdquantile(x, 0.5))

#' @title High-Density Median Absolute Deviation
#' @description Calculates the high-density (HD) median absolute deviation (MAD) of a numeric vector.
#' @param x A numeric vector.
#' @return A numeric value representing the HD MAD.
#' @examples
#' hdmad(c(1, 2, 3, 4, 5))
hdmad <- function(x) 1.4826 * hdmedian(abs(x - hdmedian(x)))

#' @title Pooled High-Density Median Absolute Deviation
#' @description Calculates the pooled HD MAD for two numeric vectors.
#' @param x A numeric vector representing the first sample.
#' @param y A numeric vector representing the second sample.
#' @return A numeric value of the pooled HD MAD.
#' @examples
#' phdmad(c(1, 2, 3), c(4, 5, 6))
phdmad <- function(x, y) pooled(x, y, hdmad)

#' @title Gamma Effect Size
#' @description Calculates the gamma effect size between two samples using high-density quantiles.
#' @param x A numeric vector representing the first sample.
#' @param y A numeric vector representing the second sample.
#' @param prob A numeric value between 0 and 1 indicating the quantile probability.
#' @return A numeric value representing the gamma effect size.
#' @examples
#' gammaEffectSize(c(1, 2, 3), c(4, 5, 6), 0.5)
gammaEffectSize <- function(x, y, prob) {
  if(length(na.exclude(y)) < 300) y <- c(na.exclude(y), rep(0, length(na.exclude(x)) - length(na.exclude(y))))
  if(length(na.exclude(x)) < 300) x <- c(na.exclude(x), rep(0, length(na.exclude(y)) - length(na.exclude(x))))
  res <- as.numeric((hdquantile(na.exclude(y), prob) - hdquantile(na.exclude(x), prob)) / phdmad(na.exclude(x), na.exclude(y)))
  return(res)
}

#' @title Count Dataset Directions
#' @description Counts the number of unique datasets with a specific `facet` and `direction`.
#' @param df Dataframe containing the data.
#' @param facet The facet variable to filter.
#' @param direction The direction to filter.
#' @return Integer value representing the count of datasets.
count_directions <- function(df, facet, direction) {
  df %>% 
    filter(facet == !!facet, direction == !!direction) %>% 
    distinct(dataset) %>% 
    nrow()
}

#' @title Calculate Effect Size
#' @description Computes the normalized effect size for a specific `facet` and `direction` 
#' using the BBGDM summary data.
#' @param bbgdm_summary Dataframe containing BBGDM results.
#' @param facet The facet variable to filter.
#' @param direction The direction to filter.
#' @return A numeric value with bootstrapped confidence intervals for the effect size.
calculate_effect_size <- function(bbgdm_summary, facet, direction) {
  result <- bbgdm_summary %>%
    filter(predictor == "hfp") %>% 
    select(facet,direction,effect_size) %>% 
    filter(facet == !!facet, direction == !!direction) %>%
    pull(effect_size) %>% 
    smean.cl.boot(conf.int=.95, B=1000, na.rm=TRUE, reps=FALSE)
  
  return(result)
}


#' @title Convert Text Size
#' @description Converts the text size in points to mm to match `geom_text` and `annotate` with theme specifications.
#' @param size Numeric value representing the text size in points.
#' @return Numeric value of size converted to millimeters.
convert_size <- function(size) {
  size <- size * 0.3528 
  return(size)
}


#' @title Create GDM Data Table
#' @description Prepares a GDM data table from dissimilarity data and predictor variables.
#' @param x Matrix or dataframe of dissimilarity values.
#' @param pred Dataframe of predictor variables.
#' @return A dataframe formatted as GDM input with class `gdmData`.
make_tab <- function(x, pred) {
  gdmTab <- x %>%
    as.matrix() %>%
    as.data.frame.table %>%
    magrittr::set_names(c("s2", "s1", "distance")) %>%
    dplyr::mutate(across(c(s2, s1), as.character)) %>%
    dplyr::left_join(pred %>%
                       magrittr::set_names(paste0("s1.", names(pred))), by = c("s1" = "s1.site")) %>%
    dplyr::left_join(pred %>%
                       magrittr::set_names(paste0("s2.", names(pred))), by =
                       c("s2" = "s2.site")) %>%
    dplyr::mutate(across(s2:s1, ~ gsub("-","_",.x))) %>% 
    dplyr::mutate(grp = paste0(pmin(s1, s2), "-", pmax(s1, s2))) %>%
    dplyr::distinct(grp, .keep_all = TRUE) %>%
    dplyr::filter(s1 != s2) %>%
    dplyr::mutate(weights = 1) %>%
    dplyr::relocate(s2, s1, grp, distance, weights, s1.x, s1.y, s2.x, s2.y) %>%
    dplyr::rename(
      s1.xCoord = "s1.x",
      s2.xCoord = "s2.x",
      s1.yCoord = "s1.y",
      s2.yCoord = "s2.y"
    ) %>%
    tibble::column_to_rownames("grp") %>%
    dplyr::select(-s2, -s1) %>%
    tidyr::drop_na()
  class(gdmTab) <- c("gdmData", "data.frame")
  return(gdmTab)
}


#' @title Run BBGDM Models
#' @description Runs Bayesian Bootstrap Generalized Dissimilarity Models (BBGDM) for differentiation and homogenization.
#' @param x Matrix or dataframe of dissimilarity values.
#' @param pred Dataframe of predictor variables.
#' @return A list containing BBGDM models for differentiation and homogenization.
run_bbgdm <- function(x, pred) {
  data_dissim <- make_tab(x, pred)
  data_sim <- data_dissim %>%  mutate(distance = 1 - distance)
  
  options(warn = 0) # turn off warnings
  model_dissim <-
    bbgdm(
      data = data_dissim ,
      geo = TRUE,
      knots = NULL,
      splines = NULL,
      boot = 1000,
      ncores = 6
    )
  model_sim <-
    bbgdm(
      data = data_sim,
      geo = TRUE,
      knots = NULL,
      splines = NULL,
      boot = 1000,
      ncores = 6
    )
  
  model <-
    list(differentiation = model_dissim, homogenization = model_sim)

 return(model)
}

#### LOAD CUSTOM BBGDM FUNCTIONS FROM SKIP
############
#' @title Running a bbgdm within the gdm modelling framework.
#' @rdname bbgdm
#' @name bbgdm
#' @description creates a \code{bbgdm} model, comprising multiple
#'   gdms. \code{bbgdm} models are core of \code{bbgdm}.
#' @param \dots for \code{bbgdm()}: one or more \code{plot()}, \code{print()} and
#'   \code{predict()}: further arguments passed to or from other methods
#' @param gdmTab formula for bbgdm model
#' @param geo presence absence matrix, sp as columns sites as rows.
#' @param splines environmental or spatial covariates at each site.
#' @param knots
#' @param bootstraps
#' @param cores
#'
#' @author Skipton Woolley
#'
#' @importFrom surveillance plapply
#' @importFrom plyr ldply
#'
#' @examples
#' library(gdm)
#' load(system.file("./data/gdm.RData", package="gdm"))
#' sppTab <- gdmExpData[, c("species", "site", "Lat", "Long")]
#' envTab <- gdmExpData[, c(2:ncol(gdmExpData))]
#' gdmTab <- formatsitepair(sppTab, bioFormat=2, XColumn="Long", YColumn="Lat",
#'                          sppColumn="species", siteColumn="site", predData=envTab)
#' fm100 <- bbgdm(gdmTab, bootstraps=100, ncores=3)
bbgdm <- function(data, geo = TRUE, splines = NULL, knots = NULL, bootstraps = 100, ncores = 1) {
    ## generate the number of sites for bayesian bootstrap.

    sitecols <- c('s1.xCoord', 's1.yCoord', 's2.xCoord', 's2.yCoord')
    sitedat <-
      rbind(as.matrix(data[, sitecols[1:2]]), as.matrix(data[, sitecols[3:4]]))
    nsites <-  data %>%
      rownames_to_column("sites") %>%
      dplyr::select(sites) %>%
      separate(sites, into = c("s1", "s2"), sep = "-") %>%
      distinct(s1) %>%
      pull() %>%
      length %>%
      {
        . + 1
      }
set.seed(123)  
    mods <-
      surveillance::plapply(
        1:bootstraps,
        bb_apply,
        nsites,
        data,
        geo,
        splines,
        knots,
        .parallel = ncores,
        .verbose = FALSE
      )
    intercept <- sapply(mods, function(x) ifelse(!is.null(x$intercept),x$intercept, 0)) %>% unlist()

    make_coef_table <- function(tab) {
      foo<-matrix(tab,ncol=3,nrow=6, byrow=TRUE,
                  dimnames=list(c("Geographic",gsub("s1.","",names(data)[7:11])),
                                paste0("coefficient.",1:3)))
      return(foo)
    }
    
    coefs<-lapply(mods, function(x) {
      if(!is.null(x$coefficients)){ data.frame(make_coef_table(tab=x$coefficients))} else NULL#{data.frame(make_coef_table(tab=0))}
    })
    
    r2 <-  sapply(mods, function(x) ifelse(!is.null(x$explained),x$explained, 0)) %>% unlist()
    
    results <- list(
     #mods=mods, #we dont want the models because they make the file too heavy
     #intercept = intercept, #we dont want the intercepts because they are not relevant to us.
      coefs = coefs,
      r2 = r2
    )
    class(results) <- c("list")
    return(results)
}

#' @title BBGDM Application Function
#' @description Creates weights for Bayesian bootstrap and fits a GDM model for each replicate.
#' @param x Bootstrap replicate number.
#' @param nsites Number of sites in the dataset.
#' @param data Dataframe containing site-pair information.
#' @param geo Logical indicating whether to include geographic coordinates as predictors.
#' @param splines Number of splines for predictors.
#' @param knots Knot locations for splines.
#' @return A GDM model object for the replicate.
bb_apply <- function(x, nsites, data, geo, splines, knots) {
  # creates the weights for bayes boot
  w <- gtools::rdirichlet(nsites, rep(1 / nsites, nsites))
  wij <- w %*% t(w)
  wij <- wij[upper.tri(wij)]
  tmp <- data
  tmp$weights <- wij
  x <- gdm::gdm(tmp,
                geo = geo,
                splines = splines,
                knots = knots)
  return(x)
}

#' @title Running a bbgdm within the gdm modelling framework.
#' @rdname bbgdm
#' @name bbgdm.predict
#' @description creates a \code{bbgdm} model, comprising multiple
#'   gdms. \code{bbgdm} models are core of \code{bbgdm}.
#' @param mod bbgdm model object
#' @param newobs new sitepair table observations, can calculate this using the gdm functions.
#' @param bbgdm is
#'
#' @author Skipton Woolley
#'
#' @importFrom surveillance plapply
#' @importFrom plyr ldply
#'
#' @examples
#' preds <- predict.bbgdm(fm100,gdmTab)
#' ## sanity check, should be pretty correlated
#' plot(gdmTab$distance,preds)
predict.bbgdm <- function(mod, newobs) {
  predicted <- rep(0, times = nrow(newobs))
  object <- mod$gdms[[1]]
  knots <- mod$gdms[[1]]$knots
  splineNo <- mod$gdms[[1]]$splines
  intercept <- mod$median.intercept
  coefs <- mod$median.coefs
  
  pred <- .C(
    "GDM_PredictFromTable",
    as.matrix(newobs),
    as.integer(FALSE),
    as.integer(length(object$predictors)),
    as.integer(nrow(newobs)),
    as.double(knots),
    as.integer(splineNo),
    as.double(c(intercept, coefs)),
    preddata = as.double(predicted),
    PACKAGE = "gdm"
  )
  
  pred$preddata
  
}

#' @rdname bbgdm
#'
#' @method plot bbgdm
#' @param add_coefs if TRUE plot coefficent estimates on prediction plot
#' @param plot_derivate if TRUE plot derivated (numberically estimated) of spline predictions
plot.bbgdm <- function(mods, add_coefs = FALSE, plot_derivate = FALSE,...) {
    gdms <- mods$gdms
    if (plot_derivate == FALSE) {
      PSAMPLE <- 200
      preddata <- rep(0, times = PSAMPLE)
      preds <- length(gdms[[1]]$predictors)
      predmax <- 0
      splineindex <- 1
      numsplines <- 3
      
      for (i in 1:preds) {
        predplotdat <- lapply(gdms, function(x)
          .C(
            "GetPredictorPlotData",
            pdata = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(x$coefficients[splineindex:(splineindex +
                                                    numsplines - 1)]),
            as.double(x$knots[splineindex:(splineindex +
                                             numsplines - 1)]),
            as.integer(numsplines),
            PACKAGE = "gdm"
          ))
        
        ## generating posterior stats from runs.
        quant.preds <-
          apply(plyr::ldply(predplotdat, function(x)
            c(x$pdata)), 2,
            function(x)
              quantile(x, c(.05, .5, .95), na.rm = T))
        varNam <- gdms[[1]]$predictors
        v <- max(quant.preds)
        if (v > predmax)
          predmax <- v
        
        splineindex <- splineindex + numsplines
      }
      
      splineindex <- 1
      numsplines <- 3
      for (i in 1:preds) {
        predplotdat <- lapply(gdms, function(x)
          .C(
            "GetPredictorPlotData",
            pdata = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(x$coefficients[splineindex:(splineindex +
                                                    numsplines - 1)]),
            as.double(x$knots[splineindex:(splineindex +
                                             numsplines - 1)]),
            as.integer(numsplines),
            PACKAGE = "gdm"
          ))
        
        ## generating posterior stats from runs.
        quant.preds <-
          apply(plyr::ldply(predplotdat, function(x)
            c(x$pdata)), 2,
            function(x)
              quantile(x, c(.05, .5, .95), na.rm = T))
        varNam <- gdms[[1]]$predictors[i]
        
        env_grad <-
          seq(
            from = gdms[[1]]$knots[splineindex],
            to = gdms[[1]]$knots[(splineindex + numsplines - 1)],
            length = PSAMPLE
          )
        plot(
          env_grad,
          quant.preds[2, ],
          xlab = varNam,
          ylab = paste("f(", varNam, ")", sep = ""),
          ylim = c(0, predmax),
          type = "n"
        )
        polygon(c(rev(env_grad), env_grad),
                c(rev(quant.preds[3, ]), quant.preds[1, ]),
                col = 'grey90',
                border = NA)
        lines(env_grad, quant.preds[2, ], col = 'dodgerblue', lwd = 2)
        
        if (add_coefs == TRUE) {
          quantcoefs <-
            apply(plyr::ldply(x, function(y)
              c(y$coefficients)), 2, function(y)
                quantile(y, c(.05, .5, .95), na.rm = T))
          
          knots <-
            gdms[[1]]$knots[splineindex:(splineindex + numsplines - 1)]
          # hack: we draw arrows but with very special "arrowheads"
          par(new = T)
          plot(
            knots,
            quantcoefs[2, splineindex:(splineindex + numsplines - 1)],
            pch = 16,
            col = 'tomato',
            axes = F,
            xlab = NA,
            ylab = NA,
            ylim = c(0, max(quantcoefs))
          )
          arrows(
            knots,
            quantcoefs[1, splineindex:(splineindex + numsplines - 1)],
            knots,
            quantcoefs[3, splineindex:(splineindex + numsplines - 1)],
            col = 'tomato',
            length = 0,
            angle = 90,
            code = 3,
            lwd = 2
          )
          axis(
            side = 4,
            col = "tomato",
            col.axis = "tomato",
            las = 1
          )
          mtext(
            "Spline coefficent estimates",
            side = 4,
            col = "tomato",
            line = -1.5,
            outer = TRUE
          )
        }
        splineindex <- splineindex + numsplines
      }
      
      
      
      
      
      
    } else {
      PSAMPLE <- 200
      preddat <- rep(0, times = PSAMPLE)
      preds <- length(gdms[[1]]$predictors)
      predmin <- 0
      predmax <- 0
      splineindex <- 1
      numsplines <- 3
      
      for (i in 1:preds) {
        predplotdat <- lapply(gdms, function(x)
          .C(
            "GetPredictorPlotData",
            pdata = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(x$coefficients[splineindex:(splineindex +
                                                    numsplines - 1)]),
            as.double(x$knots[splineindex:(splineindex +
                                             numsplines - 1)]),
            as.integer(numsplines),
            PACKAGE = "gdm"
          ))
        
        
        diffquants <-
          apply(apply(plyr::ldply(predplotdat, function(x)
            c(x$pdata)), 1, diff), 1, function(x)
              quantile(x, c(.1, .5, .9)))
        
        varNam <- gdms[[1]]$predictors[i]
        vmin <- min(diffquants)
        vmax <- max(diffquants)
        if (vmin < predmin)
          predmin <- vmin
        if (vmax > predmax)
          predmax <- vmax
        splineindex <- splineindex + numsplines
      }
      
      splineindex <- 1
      numsplines <- 3
      for (i in 1:preds) {
        predplotdat <- lapply(gdms, function(x)
          .C(
            "GetPredictorPlotData",
            pdata = as.double(preddata),
            as.integer(PSAMPLE),
            as.double(x$coefficients[splineindex:(splineindex +
                                                    numsplines - 1)]),
            as.double(x$knots[splineindex:(splineindex +
                                             numsplines - 1)]),
            as.integer(numsplines),
            PACKAGE = "gdm"
          ))
        
        
        diffquants <-
          apply(apply(plyr::ldply(predplotdat, function(x)
            c(x$pdata)), 1, diff), 1, function(x)
              quantile(x, c(.1, .5, .9)))
        
        varNam <- gdms[[1]]$predictors[i]
        
        env_grad <-
          seq(
            from = gdms[[1]]$knots[splineindex],
            to = gdms[[1]]$knots[(splineindex + numsplines - 1)],
            length = PSAMPLE
          )
        
        #calc approximate derivates.
        dX <- diff(env_grad)
        dYl <- diffquants[1, ] / dX
        dYm <- diffquants[2, ] / dX
        dYu <- diffquants[3, ] / dX
        dX <- rowMeans(embed(env_grad, 2))
        plot(
          dX,
          dYm,
          xlab = varNam,
          ylab = paste("f'(", varNam, ")", sep = ""),
          type = "n",
          ylim = c(min(dYl), max(dYu))
        )
        polygon(c(rev(dX), dX),
                c(rev(dYu), dYl),
                col = 'grey90',
                border = NA)
        lines(dX, dYm, col = 'tomato', lwd = 2)
        splineindex <- splineindex + numsplines
      }
    }
  }

#' @title Transform Environmental Data Using BBGDM
#' @description Transforms environmental data using a BBGDM model.
#' @param model A BBGDM model object.
#' @param data A raster or dataframe of environmental data.
#' @return A transformed raster stack or dataframe.
bbgdm.transform <- function(model, data) {
  #################
  ##lines used to quickly test function
  #model <- gdmModel
  #data <- climCurrExt
  #data <- cropRasts[[3:nlayers(cropRasts)]]
  #################
  options(warn.FPU = FALSE)
  rastDat <- NULL
  dataCheck <- class(data)
  
  ##error checking of inputs
  ##checks to make sure a gdm model is given
  if (class(model)[1] != "bbgdm") {
    stop("model argument must be a gdm model object")
  }
  ##checks to make sure data is a correct format
  if (!(
    dataCheck == "RasterStack" |
    dataCheck == "RasterLayer" |
    dataCheck == "RasterBrick" | dataCheck == "data.frame"
  )) {
    stop("Data to be transformed must be either a raster object or data frame")
  }
  
  ##checks rather geo was T or F in the model object
  geo <- model$gdms[[1]]$geo
  
  ##turns raster data into dataframe
  if (dataCheck == "RasterStack" |
      dataCheck == "RasterLayer" | dataCheck == "RasterBrick") {
    ##converts the raster object into a dataframe, for the gdm transformation
    rastDat <- data
    data <- rasterToPoints(rastDat)
    ##determines the cell number of the xy coordinates
    rastCells <- cellFromXY(rastDat, xy = data[, 1:2])
    
    ##checks for NA in the
    checkNAs <- as.data.frame(which(is.na(data), arr.ind = T))
    if (nrow(checkNAs) > 0) {
      warning(
        "After extracting raster data, NAs found from one or more layers. Removing NAs from data object to be transformed."
      )
      data <- na.omit(data)
      rastCells <- rastCells[-c(checkNAs$row)]
    }
    
    ##if geo was not T in the model, removes the coordinates from the data frame
    if (geo == FALSE) {
      data <- data[, 3:ncol(data)]
    }
  }
  
  sizeVal <- 10000000
  ##sets up the data to be transformed into pieces to be transformed
  holdData <- data
  fullTrans <- matrix(0, nrow(holdData), ncol(holdData))
  rows <- nrow(holdData)
  istart <- 1
  iend <- min(sizeVal, rows)
  ##to prevent errors in the transformation of the x and y values when geo is a predictor,
  ##extracts the rows with the minimum and maximum x and y values, these rows will be added
  ##onto the "chuck" given to transform, and then immediately removed after the transformation,
  ##this makes sure that the c++ code will always have access to the minimum and maximum
  ##x and y values
  if (geo == TRUE) {
    if (dataCheck == "RasterStack" |
        dataCheck == "RasterLayer" | dataCheck == "RasterBrick") {
      xMaxRow <- holdData[which.max(holdData[, "x"]), ]
      xMinRow <- holdData[which.min(holdData[, "x"]), ]
      yMaxRow <- holdData[which.max(holdData[, "y"]), ]
      yMinRow <- holdData[which.min(holdData[, "y"]), ]
    }
  }
  
  ##transform the data based on the gdm
  ##part of a loop to prevent memory errors
  while (istart < rows) {
    ##Call the dll function
    data <- holdData[istart:iend, ]
    ##adds coordinate rows to data to be transformed
    if ((dataCheck == "RasterStack" |
         dataCheck == "RasterLayer" |
         dataCheck == "RasterBrick") & geo == TRUE) {
      data <- rbind(xMaxRow, xMinRow, yMaxRow, yMinRow, data)
    }
    transformed <- matrix(0, nrow(data), ncol(data))
    z <- .C(
      "GDM_TransformFromTable",
      as.integer(nrow(data)),
      as.integer(ncol(data)),
      as.integer(model$gdms[[1]]$geo),
      as.integer(length(model$gdms[[1]]$predictors)),
      as.integer(model$gdms[[1]]$splines),
      as.double(model$gdms[[1]]$knots),
      as.double(model$median.coefs),
      as.matrix(data),
      trandata = as.double(transformed),
      PACKAGE = "gdm"
    )
    
    ## Convert transformed from a vector into a dataframe before returning...
    nRows <- nrow(data)
    nCols <- ncol(data)
    
    ## z$trandata is the transformed data vector created
    myVec <- z$trandata
    pos <- 1
    ##fills out dataframe with transformed values
    for (i in seq(from = 1, to = nCols, by = 1)) {
      tmp <- myVec[seq(from = pos, to = pos + nRows - 1)]
      transformed[, i] <- tmp
      pos <- pos + nRows
    }
    
    ##remove the coordinate rows before doing anything else
    if ((dataCheck == "RasterStack" |
         dataCheck == "RasterLayer" |
         dataCheck == "RasterBrick") & geo == TRUE) {
      transformed <- transformed[-c(1:4), ]
    }
    
    ##places the transformed values into the readied data frame
    fullTrans[istart:iend, ] <- transformed
    istart <- iend + 1
    iend <- min(istart + (sizeVal - 1), rows)
  }
  
  ##if wanted output data as raster, provides maps raster, or output table
  if (dataCheck == "RasterStack" |
      dataCheck == "RasterLayer" | dataCheck == "RasterBrick") {
    ##maps the transformed data back to the input rasters
    rastLay <- rastDat[[1]]
    rastLay[] <- NA
    outputRasts <- stack()
    for (nn in 1:ncol(fullTrans)) {
      #print(nn)
      #nn=1
      holdLay <- rastLay
      holdLay[rastCells] <- fullTrans[, nn]
      #holdLay[rastCells] <- holdData[,nn]
      
      outputRasts <- stack(outputRasts, holdLay)
    }
    ##renames raster layers to be the same as the input
    if (geo) {
      names(outputRasts) <- c("xCoord", "yCoord", names(rastDat))
    } else {
      names(outputRasts) <- names(rastDat)
    }
    
    ##get the predictors with non-zero sum of coefficients
    splineindex <- 1
    predInd <- NULL
    for (i in 1:length(model$gdms[[1]]$predictors)) {
      #i <- 1
      ##only if the sum of the coefficients associated with this predictor is > 0.....
      numsplines <- model$gdms[[1]]$splines[i]
      if (sum(model$median.coefs[splineindex:(splineindex + numsplines -
                                              1)]) > 0) {
        predInd <- c(predInd, i)
      }
      splineindex <- splineindex + numsplines
    }
    if (geo) {
      predInd <- c(1, 2, predInd[-1] + 1)
    }
    
    outputRasts <- outputRasts[[predInd]]
    
    ##returns rasters
    return(outputRasts)
  } else{
    if (is.null(rastDat)) {
      ##if not raster data, sends back the transformed data
      colnames(fullTrans) <- colnames(data)
      return(fullTrans)
    } else{
      ##returns only the transformed variable data as a table, and the cells with which to map to
      colnames(fullTrans) <- colnames(data)
      return(list(fullTrans, rastCells))
    }
  }
}

#' @title Compute BBGDM Residuals
#' @description Computes random quantile residuals for a BBGDM model.
#' @param object A BBGDM model object.
#' @return A list containing residuals, predicted probabilities, and observed dissimilarities.
residuals.bbgdm <- function (object) {
  y <- object$gdms[[1]]$observed  # offset <- object$offset
  preds <-
    t(plyr::ldply(object$gdms, function(x) {
      x$predicted
    }, .progress = "none"))
  pi <- apply(preds, 1, mean)
  a <- pbinom(y - 1, 1, pi)#-1
  b <- pbinom(y, 1, pi)
  u <- runif(n = length(y), min = a, max = b)
  # u[u==1] <- u[u==1]-1e-5
  # u[u==0] <- u[u==0]+1e-5
  res <- qnorm(u)
  structure(list(res = res, pi = pi, y = y))
}

#' @title Compute BBGDM Residuals
#' @description Computes random quantile residuals for a BBGDM model.
#' @param object A BBGDM model object.
#' @return A list containing residuals, predicted probabilities, and observed dissimilarities.
plot.residuals <- function(x, ...) {
  qqnorm(x$res, ylab = 'Random Quantile Residuals', main = "")
  qqline(x$res, col = 'red')
  hist(x$res, xlab = "Random Quantile Residuals", main = "", ...)
  plot(x$res, x$pi, xlab = "Predicted Dissimilarity", ylab = "Random Quantile Residuals", ...)
  plot(x$pi, x$y, xlab = "Predicted Dissimilarity", ylab = "Observed Dissimilarity", ...)
}

#' @title BBGDM Wald Test
#' @description Performs a Wald test for coefficients in a BBGDM model.
#' @param object A BBGDM model object.
#' @param H0 Null hypothesis value (default is 0).
#' @return A matrix of Wald statistics and p-values for model coefficients.
bbgdm.wald.test <- function(object, H0 = 0) {
  # H0: Hypothesis test = 0
  # IM: Identiy Matrix
  # beta: parameter estimates from model
  # vcov: Variance-covariance matrix estimated from BB
  index <-
    which(!sapply(object$gdms, is.null, simplify = TRUE))[1] #retrieve index for baseline model
  A <-
    plyr::ldply(object$gdms, function(x) {
      c(x$intercept, x$coefficients)
    }, .progress = "none")# all.coefs.se #matrix of B bootstrap coeficient estimates.
  esti.var <- var(A) #make sure this is a matrix
  beta <-
    c(object$median.intercept, object$median.coefs) #medians of coef estimates
  
  
  
  #Intercept
  intercept_IM <- matrix(c(1, rep(0, length(beta) - 1)), nrow = 1)
  wd_inter <-
    t(intercept_IM %*% beta - H0) %*% solve(intercept_IM %*% esti.var %*% t(intercept_IM)) %*%
    (intercept_IM %*% beta - H0)
  pval_i = 1 - pchisq(wd_inter, 1)
  
  #Splines
  # beta <- c(object$median.coefs) #medians of coef estimates
  splineLength <- object$gdms[[index]]$splines
  val1 <- seq(2, length(beta), splineLength[1])
  val2 <- seq(1 + splineLength[1], length(beta), splineLength[1])
  wd_vals <- matrix(NA, length(val1) + 1, 3)
  wd_vals[1, ] <- c(wd_inter, 1, pval_i)
  for (i in 1:length(splineLength)) {
    w <- splineLength[1]
    L <- matrix(rep(0, length(beta) * w), ncol = length(beta))
    Terms <- seq(val1[i], val2[i], 1)
    for (ii in 1:w)
      L[ii, Terms[ii]] <- 1
    vcov2 <- try(solve(L %*% esti.var %*% t(L)), silent = TRUE)
    if ((class(vcov2) == "try-error") || (class(vcov2) == "try-error"))
    {
      wd <- 0
    } else{
      wd <- t(L %*% beta - H0) %*% vcov2 %*% (L %*% beta - H0)
    }
    pv <- 1 - pchisq(wd, df = w)
    wd_vals[1 + i, ] <- c(wd, w, pv)
  }
  colnames(wd_vals) <- c("bbgdm_W", "bbgdm_df", "bbgdm_p-value")
  rownames(wd_vals) <- c('intercept', object$gdms[[index]]$predictors)
  wd_vals
}

#' @title Negative Exponential Link Function
#' @description Provides a negative exponential link function for GLMs.
#' @return A list representing the negative exponential link function.
negexp <- function() {
  linkfun <- function(mu)
    - log(1 - mu)
  linkinv <- function(eta)
    1 - exp(-eta)
  mu.eta <- function(eta)
    exp(-eta)
  valideta <- function(eta)
    all(is.finite(eta))
  link <- paste0("negexp")
  structure(
    list(
      linkfun = linkfun,
      linkinv = linkinv,
      mu.eta = mu.eta,
      valideta = valideta,
      name = link
    ),
    class = "link-glm"
  )
}