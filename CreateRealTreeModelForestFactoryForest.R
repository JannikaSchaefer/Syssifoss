#### load libraries and data ####
library(data.table)
library(lidR)

## Forest Factory output
forestInfo <- fread("2022-11-22_13-29-26_471286_fertile_soil_8pft.res", skip = 2)


## single tree point clouds
pcPath <- "./pointclouds/"
treeModelFiles <- list.files(pcPath, pattern = ".laz")


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

#### Data modification ####
### add species information
species <- data.frame(id = 1:8, Spec = c("PinSyl", "PicAbi", "FagSyl", "Quercus", "SorAuc", "PopTre", "BetPub", "BetPen")) 

forestInfo <- merge(forestInfo, species, by.x = "Grp", by.y = "id")

### add unique ID for Plots
forestInfo$PlotHec <- paste0("Plot", forestInfo$Plot, "_Hec", forestInfo$Hec)

### keep only interesting columns
forests <- forestInfo[, c("Grp", "Time", "Species", "BT", "D", "H", "AGE" , "Plot", "Hec", "X", "Y", "CLP", "LAI", "CD", "ID", "CR", "Spec", "PlotHec")]

### shift tree coordinates to range from 0 to 100
forests[, Xorig := X]
forests[, Yorig := Y]
forests[, X := Xorig - (Xorig %/% 100) * 100]
forests[, Y := Yorig - (Yorig %/% 100) * 100]


#### Load tree model data ####
treeModelInfoMetrics <- fread("./csv/result_metrics.csv")
treeModelInfoGeneral <- fread("./csv/result_general.csv" )
treeModelInfo <- merge(treeModelInfoMetrics, treeModelInfoGeneral)

## extract species short name
treeModelInfo$Species <- substr(treeModelInfo$tree_id, 1, 6)

## change "QuePet" and "QueRob" to "Quercus"
treeModelInfo$Species[treeModelInfo$Species %in% c("QuePet", "QueRob")] <- "Quercus"

###filter tree models according to tree model requirements 
## define tree model requirements
sources <- "ULS" 
canopy <- "leaf-on"
mtType <- "ULS-on"
# qualities <- ("q1", "q2", "q3", "q4", "q5", "q6")[1:3] ## quality information is not yet included in the csv file

treeModels <- treeModelInfo[treeModelInfo$source %in% sources & treeModelInfo$canopy_condition %in% canopy, ] #treeModelInfo[treeModelInfo$source %in% sources & treeModelInfo$canopy_condition %in% canopy & treeModelInfo$quality %in% qualities, ]


#### Select tree models ####
hectars <- unique(forests$Hec)
for(hectar in hectars){
  stand <- forests[Hec == hectar]
  for(i1 in 1:nrow(stand)){
    tree <- stand[i1, ]
    ## select for matching species
    treeModels1 <- treeModels[Species == tree$Spec]
    ## select for matching height
    treeModels2 <- treeModels1[height_m >= tree$H - 4 & height_m <= tree$H + 4]
    if(nrow(treeModels2) == 0){
      selTree <- treeModels1[which(abs(tree$H - treeModels1$height_m) == min(abs(tree$H - treeModels1$height_m)))]
    } else {
      selTree <- treeModels2[sample(1:nrow(treeModels2), 1), ]
    }
    
    stand[i1, c("mtTreename", "tmSpecies", "mtH", "mtCBH", "mtCD", "mtPositionX", "mtPositionY", "mtPositionZ", "canopyCondition", "tmSource", "tmDate")] <-  selTree[1, .(tree_id, Species, height_m, crown_base_height_m, mean_crown_diameter_m, x, y, z, canopy_condition, source, date)]
    
  }
  
  #### Create synthetic forest stand ####
  stand$TreeId <- 1:nrow(stand) 
  
  ## calculate scaling factor for height
  stand$z_scale <- stand$H/stand$mtH
  
  ## random rotation angle
  stand$random_rot <- sample(seq(1,360,1), nrow(stand), replace = T)
  
  
  #### Create a tree point cloud forest stand ####
  forestList <- list()
  
  for(i1 in 1:nrow(stand)){
    ## load tree point cloud
    treename <- paste0(stand$mtTreename[i1], "_.*", mtType, ".laz")
    treeModelFileName <- treeModelFiles[grep(treename, treeModelFiles)]
    modeltreeLAS <- readLAS(paste0(pcPath, treeModelFileName), select = "xyz")
    tm <- as.data.frame(modeltreeLAS@data)
    
    ## shift tree point cloud position to 0,0,0
    tm$shiftX <- tm$X - stand$mtPositionX[i1] 
    tm$shiftY <- tm$Y - stand$mtPositionY[i1] 
    tm$shiftZ <- tm$Z - stand$mtPositionZ[i1]
    
    ## scale in z direction -> height
    tm$scaledZ <- tm$shiftZ * stand$z_scale[i1]
    
    ## Scale in xy direction 
    tm$scaledX <- tm$shiftX * stand$z_scale[i1]
    tm$scaledY <- tm$shiftY * stand$z_scale[i1]
    
    
    
    ## rotate tree point cloud around z axis
    rotXY <- turn_mat(stand$random_rot[i1], tm$scaledX, tm$scaledY)
    tm$rotX <- rotXY[, 1]
    tm$rotY <- rotXY[, 2]
    
    ## shift tree point cloud position to Forest Factory tree position
    tm$finalX <- tm$rotX + stand$X[i1]
    tm$finalY <- tm$rotY + stand$Y[i1] 
    tm$finalZ <- tm$scaledZ + 0
    
    
    ## add additional information
    tm$id <- stand$TreeId[i1]
    tm$Spec <- stand$Spec[i1]
    tm$SpecID <- stand$Species[i1]
    tm$Plot <- stand$Plot[i1]
    tm$Hec <- stand$Hec[i1]
    tm$PlotHec <- stand$PlotHec[i1]
    
    ## Save point cloud
    forestList[[i1]] <- tm
    
   
  }
  
  #### Export tree point cloud forest as txt file ####
  tpcForest <- rbindlist(forestList)
  out <- tpcForest[, c("finalX", "finalY", "finalZ", "id", "Spec", "SpecID", "Hec", "PlotHec")]
  
  outName1 <- paste0("PC_FF_", hectar, ".txt")
  writeLines("//\"X\" \"Y\" \"Z\" \"id\" \"Spec\" \"SpecID\" \"Plot\" \"Hec\" \"PlotHec\"", con = outName1)
  fwrite(out, outName1, append = T, col.names =  F, row.names = F, sep = " ")
  
  outName2 <- paste0("Info_FF_", hectar, ".txt")
  fwrite(stand, outName2)
  
}
