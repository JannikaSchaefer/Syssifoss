#### load libraries and data ####
library(data.table)
library(lidR)

## selected tree models 
forest <- fread("TreeFilteringResults.txt")

## single tree point cloud path
pcPath <- "./pointclouds/"
treeModelFiles <- list.files(pcPath, pattern = ".laz")

## define output name
outName <- "RTM_Forest.txt"

#### function to rotate point clouds around the z axis ####
turn_mat <- function(al, x, y) {
  xcen <- 0
  ycen <- 0
  
  x <- x - xcen
  y <- y - ycen
  
  newx <- cos(al*pi/180) * x - sin(al*pi/180) * y
  newy <- sin(al*pi/180) * x + cos(al*pi/180) * y
  
  newx2 <- newx + xcen
  newy2 <- newy + ycen
  
  return(cbind(newx2, newy2))
}

#### prepare data ####

## calculate scaling factor for height
forest$z_scale <- forest$Height/forest$mtH

## calculate scaling factor for crown diameter
forest$xy_scale <- forest$meanCD/forest$mtCD

## random rotation angle
forest$random_rot <- sample(seq(1,360,1), nrow(forest), replace = T)


#### insert tree model point clouds to the forest ####
forestList <- list()

for(t1 in 1:nrow(forest)){
  tree <- forest[t1, ]
  
  ## load tree point cloud
  canopyCondition <- strsplit(tree$canopyCondition, "-")[[1]][2]
  treeName <- paste0(tree$mtTreename, "_", tree$mtDate, ".*_", tree$mtSource, "-", canopyCondition, ".laz")
  treeModelFileName <- treeModelFiles[grep(treeName, treeModelFiles)]
  
  treeModelLAS <- readLAS(paste0(pcPath, treeModelFileName))
  tm <- as.data.frame(treeModelLAS@data)
  
  ## shift tree point cloud position to 0,0,0
  tm$shiftX <- tm$X - tree$mtPositionX 
  tm$shiftY <- tm$Y - tree$mtPositionY 
  tm$shiftZ <- tm$Z - tree$mtPositionZ 
  
  ## scale in z direction -> height
  tm$scaledZ <- tm$shiftZ * tree$z_scale
  
  ## scale in xy direction -> crown diameter
  tm$scaledX <- tm$shiftX * tree$xy_scale
  tm$scaledY <- tm$shiftY * tree$xy_scale
  
  ## rotate tree point cloud around z axis
  rotXY <- turn_mat(tree$random_rot, tm$scaledX, tm$scaledY)
  tm$rotX <- rotXY[, 1]
  tm$rotY <- rotXY[, 2]
  
  ## shift tree point cloud position to real tree position
  tm$finalX <- tm$rotX + tree$realX 
  tm$finalY <- tm$rotY + tree$realY 
  tm$finalZ <- tm$scaledZ + 0
  
  ## add id
  tm$id <- t1
  
  ## add species
  tm$species <- tree$Spec
  
  ## add tree point cloud to forest
  forestList[[t1]] <- tm[, c("id", "finalX", "finalY", "finalZ", "species")]
  rm(tree, canopyCondition, treeName, treeModelFileName, treeModelLAS, tm, rotXY)
}


#### Export tree point cloud forest as txt file ####
tpcForest <- rbindlist(forestList)
writeLines("//\"X\" \"Y\" \"Z\" \"id\"", con = outName)
fwrite(tpcForest, outName, append = T, col.names =  F, row.names = F, sep = " ")
