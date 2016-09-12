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

