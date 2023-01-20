#### load libraries and data ####
library(data.table)

## selected tree models 
forest <- fread("TreeFilteringResults.txt")

## define output name
outName <- "STM_Forest.txt"

## define planar point density of the crown and maximum distance between stem points
pointDensity <- 1500 # pts/m2
maxDist <- 0.04 # in m

#### insert simplified tree point clouds to the forest ####
forestList <- list()

for(t1 in 1:nrow(forest)){
  tree <- forest[t1, ]
  
  cd <- tree$meanCD ## in m
  th <- tree$Height ## in m
  cbh <- tree$CBH2 ## crown base height in m
  dbh <- tree$DBH * 0.01 ## convert cm to m
  
  ## trees without crown -> only stem
  if(is.na(cd) | cd == 0){
    ## creates a cylindrical point cloud with a diameter = tree dbh, height = tree height, and a maximum distance between points of maxDist (0.04 m)
    xstvals <- seq(-dbh/2, dbh/2, by = maxDist)
    ystvals1 <- sqrt((dbh/2)^2 - xstvals^2)
    ystvals2 <- -1 * ystvals1
    zstvals <- seq(0, th, maxDist)
    xystvals <- cbind(c(xstvals, xstvals), c(ystvals1, ystvals2))
    xystvals <- as.data.frame(xystvals)
    
    st <- expand.grid(zstvals, c(ystvals1, ystvals2))
    names(xystvals) <- c("X", "Y")
    names(st) <- c("Z", "Y")
    merged <- merge(xystvals, st)
    merged <- round(merged, 3)
    
    zyl <- merged[, c("X", "Y", "Z")]
    tm <- rbind(zyl, c(0, 0, 0))
    
    ## shift to real tree position
    tm$finalX <- tr$X + tree$realX
    tm$finalY <- tr$Y + tree$realY
    tm$finalZ <- tr$Z + 0
    
    ## add id
    tm$id <- t1
    
    ## add species
    tm$species <- tree$Spec
    
    ## add tree point cloud to forest
    forestList[[t1]] <- tm[, c("id", "finalX", "finalY", "finalZ", "species")]
    
  } else {
  ## trees with crown
    
    ## create crown point cloud
    bboxArea <- cd * cd
    pointNumber <- pointDensity * bboxArea
    randomPointsX <- sample(seq(-cd/2, cd/2, 0.01), pointNumber, replace = T)
    randomPointsY <- sample(seq(-cd/2, cd/2, 0.01), pointNumber, replace = T)
    randomPointsZ <- sample(seq(0, th, 0.01), pointNumber, replace = T)
    randomPoints <- data.frame(randomPointsX, randomPointsY, randomPointsZ)
    colnames(randomPoints) <- c("X", "Y", "Z")
    
    
    if(is.na(cbh)){
      ## if no crown base height is given, crown height is assumed to be 40 % of tree height
      cbh <- th - (0.4 * th)
    }
    
    ## filter for points wihtin crown ellipsoid
    ## ellipsoid equation: (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
    a <- 0.5 * cd
    c <- 0.5 * (th - cbh)
    ell <- randomPoints[which ((randomPoints$X/a)^2 + (randomPoints$Y/a)^2 + ((randomPoints$Z - (cbh + c))/c)^2 <= 1), ]
      
    
    ## create points of stem cylinder
    xstvals <- seq(-dbh/2, dbh/2, by = 0.02)
    ystvals1 <- sqrt((dbh/2)^2 - xstvals^2)
    ystvals2 <- -1 * ystvals1
    zstvals <- seq(0, cbh, 0.02)
    xystvals <- cbind(c(xstvals, xstvals), c(ystvals1, ystvals2))
    xystvals <- as.data.frame(xystvals)
    
    st <- expand.grid(zstvals, c(ystvals1, ystvals2))
    names(xystvals) <- c("X", "Y")
    names(st) <- c("Z", "Y")
    merged <- merge(xystvals, st)
    merged <- round(merged, 3)
    
    
    zyl <- merged[, c("X", "Y", "Z")]
    tm <- rbind(ell, zyl)
    
    ## shift to real tree position
    tm$finalX <- tm$X + tree$realX
    tm$finalY <- tm$Y + tree$realY
    tm$finalZ <- tm$Z + 0
    
    ## add id
    tm$id <- t1
    
    ## add species
    tm$species <- tree$Spec
    
    ## add tree point cloud to forest
    forestList[[t1]] <- tm[, c("id", "finalX", "finalY", "finalZ", "species")]
  }
}


#### Export tree point cloud forest as txt file ####
tpcForest <- rbindlist(forestList)
writeLines("//\"X\" \"Y\" \"Z\" \"id\"", con = outName)
fwrite(tpcForest, outName, append = T, col.names =  F, row.names = F, sep = " ")
