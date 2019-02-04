# macarthur work

--------------------------------------------------

## ./analysis_data/
Contains the input data for the analyses.  Script "analysis.Rnw" will automatically download the relevant files into this directory.

## ./shps/

Project contains 3 regions:

South East Asia (sea)
- Cambodia
- Laos
- Myanmar
- Thailand
- Vietnam  

Africa (africa)
- Burundi
- Democratic Republic of the Congo
- Kenya
- Mozambique
- Rwanda
- Tanzania
- Uganda
- Zambia  

South America (sa)
- Bolivia
- Colombia
- Ecuador
- Peru
- Venezuela  


1. Download global GADM files (http://biogeo.ucdavis.edu/data/gadm2.8/gadm28_levels.shp.zip)
For each region:
2. Using QGIS and the GADM global ADM2 shapefile, create a new shapefile which contains the ADM2 features for the relevant countries [eg: ./shps/sea/sea_adm2.shp]
3. Generate a grid (Vector>Research Tools>Vector Grid) by setting the extent to the new adm2 layer (select layer then click "update extents from layer"), setting resolution to 0.05 (lock 1:1 ratio) and outputting the grid as polygons [eg: ./shps/sea/sea_full_grid.shp]
4. Remove extraneous grid cells (Vector>Data Management Tools>Join attributes by location) by setting target layer to your grid, join layer to the regional adm2 layer, and ensure the "keep only matching records" option is selected [eg: ./shps/sea/sea_grid.shp]
5. Save resulting "clipped" grid layer as a CSV (right click layer, save as, then set format to Comma Separated Value) [eg: analysis_data/grid_info/sea_data.csv]
6. Remove all fields except "ID" from clipped grid (install QGIS "Table Manager" plugin if not already installed then Vector>Table Manager>Table Manager) by selected all fields except "ID", click the "Delete" and then the "Save" button
7. Convert "clipped" grid layer to points using the polygon centroid tool (Vector>Geometry Tools>Polygon centroids) by selecting the grid layer [eg: ./shps/sea/sea_points.shp]

_Note: if you do not use the file structue/names from descriptions, other scripts may not work properly._


## ./hansen/

dependencies
- linux packages: gdal/ogr
- python packages: rasterstats, pandas

1. download and mosaic hansen data  
`bash ./hansen/get_hansen_for_bbox`
2. run extracts  
`python ./hansen/run_hansen_extract.py`


## ./sciclone_data/

1. run extract job(s) on sciclone
2. run extract merge script  
`python ./sciclone_data/merge_data.py`

## /modelData

Correlograms
- graphs that plot how connected historical NDVI values are in cells at different distances from each other (e.g. how related are historical NDVI values in cells that are 5 km away from each other vs. 100km away from each other) and gives a value to represent that. 
- We use the correlograms to construct the continuous treatment variable, such that the "coefficient of relatedness" (the value on the y-axis at a given distance between any two cells shown on the x axis) is multiplied by the distance from a cell to a project location, to calculate a treatment value. This value grows over time as more projects become active. We exclude any projects as being too far away if they exceed the distance determined by when the correlogram line crosses the x-axis (and thus the "coefficient of relatedness" is equal to 0).
- Correlogram graphs are calculated using all cells in a country. They only need to be recalculated if the historical NDVI values change.

Datasets
- includes the data used to run the analytical models for MacArthur, if we have them (some were lost on the desktop where this data was stored)
- cambodia_panel_data_add.csv is the panel dataset used to create the Cambodia Infrastructure Regression Results (Table 3) in the original MacArthur working paper. **Note that the DMSP values extracted for this dataset are factually incorrect! ** The dataset with updated DMSP ntl values can be found in the MacArthur Box Sync folder, /modelData/cambodia_infra_add_oct2017.csv" but the regression results will not match any models from the Working Paper that include ntl baseline or pre-trend values.

## Scripts to Merge Spatial Data and Produce Panel Datasets

MacArthur_Paper_Analysis_CAMBODIA
- original script to build the Cambodia panel dataset
- Does not work! R functions and data sources have changed
- Replaced with script "MacArthur_Paper_Analysis_Cambodia_Thresh10_Recreate"

MacArthur_Paper_Analysis_Cambodia_Thresh10_Recreate
- replaces original script used to create panel dataset for analysis
- updated file paths when input data was moved to Box Sync folder called "MacArthur" (though some input data still exists in Git Repo as well)
- updated functions for things like TimeRangeTrend that referenced out of date R package built by Dan
- includes DMSP ntl data updated in fall 2017 (though still have old DMSP data extracts)
- never used to produce any reports, working papers, etc.


## Scripts for Model Analysis

Models_Cambodia.R 
- model analysis for Cambodia infrastructure sector
- includes many models, but Stargazer code at end was used to output results tables in the Working Paper (which also identifies which of the models we used for the working paper)






