
import os

import pandas as pd
import rasterstats as rs

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


rlist = ["sea", "sa", "africa"]
#rlist = ["sea"]

for rid in rlist:
    
    print rid
    
    vector1 = base_dir + "/shps/"+rid+"/"+rid+"_grid.shp"
    vector2 = base_dir + "/shps/"+rid+"/"+rid+"_points.shp"
    raster1 = base_dir + "/hansen/data/treecover2000/"+rid+"_mosaic.tif"
    raster2 = base_dir + "/hansen/data/loss/"+rid+"_mosaic.tif"

    # raster3 = base_dir + "/hansen/data/ltdr_ndvi_yearly_max_2000.tif"

    cats = {0: '0', 1: '1'}


    treecover2000 = rs.zonal_stats(vector1, raster1, stats='mean', geojson_out=True) 

    loss = rs.zonal_stats(vector1, raster2, categorical=True, category_map=cats, geojson_out=True)

    # ltdr = rs.zonal_stats(vector2, raster3, stats='mean', geojson_out=True) 


    tmp_data = [i['properties'] for i in treecover2000]

    percent_loss_data = []
    for i in range(len(loss)):
        tmp_props = loss[i]['properties']
        if not '1' in tmp_props:
            tmp_props['1'] = 0
        elif not '0' in tmp_props:
            tmp_props['0'] = 0
    
        tmp_percent_loss = 100 * tmp_props['1'] / (tmp_props['1'] + tmp_props['0'])
        percent_loss_data.append(tmp_percent_loss)


    # ltdr_data = [i['properties']['mean'] for i in ltdr]


    tmp_df = pd.DataFrame(tmp_data)
    tmp_df.rename(columns = {'mean':'tc00_e'}, inplace=True)

    tmp_df['tc00_e'] = tmp_df['tc00_e'].astype(int)
    tmp_df['per_loss'] = percent_loss_data
    # tmp_df['lnyx_2000e'] = ltdr_data

    #tmp_df['ad_extract'].fillna('NA', inplace=True)
    tmp_df.to_csv(base_dir + "/analysis_data/hansen_extracts/"+rid+"_extract.csv", sep=",", encoding="utf-8", index=False)

