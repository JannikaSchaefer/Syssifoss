#### load libraries and data ####
library(data.table)

## forest stand information on single tree level
forestInfo <- fread("ForestInventoryData.csv", dec = ",")

## tree model data files
treeModelInfoMetrics <- fread("./csv/result_metrics.csv")
treeModelInfoGeneral <- fread("./csv/result_general.csv" )

## define output name 
outName <- "TreeFilteringResults.txt"

## define tree model requirements and 
sources <- "ULS" 
canopy <- "leaf-on"
# qualities <- ("q1", "q2", "q3", "q4", "q5", "q6")[1:3] ## quality information is not yet included in the csv file

#### Prepare data ####
## merge tree model information
treeModelInfo <- merge(treeModelInfoMetrics, treeModelInfoGeneral)

treeModelInfo$date <- as.character(treeModelInfo$date)

### add species information of tree models
## extract species short name
treeModelInfo$Species <- substr(treeModelInfo$tree_id, 1, 6)

## define species groups
broad_leaved <- c("AceCam", "AcePse", "BetPen", "CarBet", "FagSyl", "FraExc", "JugReg", "PruAvi", "PruSer", "QuePet", "QueRob", "QueRub", "RobPse", "SalCap", "SorTor", "TilSpe")
coniferous <- c("AbiAlb", "LarDec", "PicAbi", "PinSyl", "PseMen", "TsuHet")

## add species group information
treeModelInfo$SpecGroup <- NA
treeModelInfo$SpecGroup[treeModelInfo$Species %in% broad_leaved] <- "Broad-leaved"
treeModelInfo$SpecGroup[treeModelInfo$Species %in% coniferous] <- "Coniferous"


#### filter tree models according to tree model requirements ####
treeModels <- treeModelInfo[treeModelInfo$source %in% sources & treeModelInfo$canopy_condition %in% canopy, ] #treeModelInfo[treeModelInfo$source %in% sources & treeModelInfo$canopy_condition %in% canopy & treeModelInfo$quality %in% qualities, ]


#### tree model selection ####
## prepare output file 
forest <- forestInfo
forest$mtTreename <- "NA"
forest$mtSpecies <- "NA"
forest$mtH <- 0
forest$mtCBH <- 0
forest$mtCD <- 0
forest$mtPositionX <- 0
forest$mtPositionY <- 0
forest$mtPositionZ <- 0
forest$canopyCondition <- "NA"
forest$mtSource <- "NA"
forest$mtDate <- "NA"

## run through trees
for(t1 in 1:nrow(forest)){
  tree <- forest[t1, ]
  
  ## include only tree models of matching species
  models <- treeModels[treeModels$Species == tree$Spec, ]
  
  ## if there is no tree model of matching species, include tree models of matching species group
  if(nrow(models) < 1){
    models <- treeModels[treeModels$SpecGroup == tree$SpecGroup, ]
  }

  ## calculate the Euclidean distance in terms of tree height and crown diameter for all remaining tree models with respect to the query tree
  models$hcdDist <- sqrt((models$height_m - tree$Height)^2 + (models$mean_crown_diameter_m - tree$meanCD)^2)
  
  ## randomly select a tree model with a Euclidean distance < 4 m or select the tree model with the smalles distance
  if(sum(models$hcdDist < 4) > 0){
    selTreeModel <- models[sample(which(models$hcdDist < 4), 1), ]  
  } else {
    selTreeModel <- models[which(models$hcdDist == min(models$hcdDist)), ]
  }
  
  ## save selected tree model information
  forest[t1, c("mtTreename", "mtSpecies", "mtH", "mtCBH", "mtCD", "mtPositionX", "mtPositionY", "mtPositionZ", "canopyCondition", "mtSource", "mtDate")] <-  selTreeModel[, .(tree_id, Species, height_m, crown_base_height_m, mean_crown_diameter_m, x, y, z, canopy_condition, source, date)]
}

#### Export tree model filtering results ####
fwrite(forest, outName)
