# Syssifoss - Synthetic structural remote sensing data for improved forest inventory models

Here, we provide scripts for creating synthetic 3D representations of forests using single tree point clouds. 

Input data should be forest inventory data on a single tree level including tree positions, species, height, crown diameter, crown base height and diameter at breast height (DBH). One example input file is provided (`ForestInventoryData.csv`).

Single tree point clouds can be real tree point clouds, i.e., point clouds of single trees extracted from real laser scanning data, or simplified tree models, i.e., point clouds in form of cylindrical stems and ellipsoidal crowns.

Real single tree point clouds can be downloaded from pytreedb (https://pytreedb.geog.uni-heidelberg.de/). You can preselect the source of the tree point clouds (ALS, ULS, TLS), the canopy condition (leaf-on, leaf-off) and the quality of the tree point cloud segmentation (q1-q6). 

Save the selected point clouds and export the metadata of the selected trees to csv files. You will get two folders, "csv" contains the files "result_general.csv" and "result_metrics.csv", and "pointclouds" contains all .laz-files of the trees. Even if you preselect a certain laser scanning source, you will get all point clouds of the trees for which the desired source is available.

`SelectTreeModels.R` selects a matching real tree point cloud for every tree in the forest. 

`CreateRealTreeModelForest.R` creates the forest point cloud based on the single tree point clouds selected in `SelectTreeModels.R`. 

`CreateSimplifiedTreeModelForest.R` creates a forest point cloud composed of simplified tree point clouds. For creating the simplified tree model forest, it is not necessary to download the real treepoint clouds.

The created 3D representations can be used as input for laser scanning simulations using HELIOS++ (https://github.com/3dgeo-heidelberg/helios). Example scripts can be downloaded here: https://doi.org/10.5445/IR/1000147797. 


