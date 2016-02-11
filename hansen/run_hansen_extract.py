
import os
import math

import pandas as pd
import rasterstats as rs

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


rlist = ["sea", "sa", "africa"]
#rlist = ["sea"]

for rid in rlist:
    
    print rid
    
    grid_vector = base_dir + "/shps/"+rid+"/"+rid+"_grid.shp"
    treecover2000_raster = base_dir + "/hansen/data/treecover2000/"+rid+"_mosaic.tif"
    loss_raster = base_dir + "/hansen/data/loss/"+rid+"_mosaic.tif"
    gain_raster = base_dir + "/hansen/data/gain/"+rid+"_mosaic.tif"
    lossyear_raster = base_dir + "/hansen/data/lossyear/"+rid+"_mosaic.tif"


    binary_cats = {0: '0', 1: '1'}

    year_cats = {}
    for i in range(15):
        if i == 0:
            year_cats[0] = 'loss_none'
        elif len(str(i)) == 1:
            year_cats[i] = 'loss_200' + str(i)
        else:
            year_cats[i] = 'loss_20' + str(i)


    treecover2000 = rs.zonal_stats(grid_vector, treecover2000_raster, stats='mean', geojson_out=True) 

    loss = rs.zonal_stats(grid_vector, loss_raster, categorical=True, category_map=binary_cats, geojson_out=True)

    gain = rs.zonal_stats(grid_vector, gain_raster, categorical=True, category_map=binary_cats, geojson_out=True)
   
    lossyear = rs.zonal_stats(grid_vector, lossyear_raster, categorical=True, category_map=year_cats, geojson_out=True)



    # points_vector = base_dir + "/shps/"+rid+"/"+rid+"_points.shp"
    # ltdr_raster = base_dir + "/hansen/data/ltdr_ndvi_yearly_max_2000.tif"
    # ltdr = rs.zonal_stats(vector2, raster3, stats='mean', geojson_out=True) 

    def roundxy(val, dec=None):
        if dec == None:
            dec = 2
        return math.floor(10**int(dec)*val)/(10**int(dec))


    tmp_data = [i['properties'] for i in treecover2000]

    percent_loss_data = []
    for i in range(len(loss)):
        tmp_props = loss[i]['properties']
        if not '1' in tmp_props:
            tmp_props['1'] = 0
        elif not '0' in tmp_props:
            tmp_props['0'] = 0
    
        tmp_percent_loss = 100 * tmp_props['1'] / float(tmp_props['1'] + tmp_props['0'])
        percent_loss_data.append(roundxy(tmp_percent_loss))


    percent_gain_data = []
    for i in range(len(gain)):
        tmp_props = gain[i]['properties']
        if not '1' in tmp_props:
            tmp_props['1'] = 0
        elif not '0' in tmp_props:
            tmp_props['0'] = 0
    
        tmp_percent_gain = 100 * tmp_props['1'] / float(tmp_props['1'] + tmp_props['0'])
        percent_gain_data.append(roundxy(tmp_percent_gain))



    percent_lossyear_data = [[] for i in range(1,15)]
    lossyear_yearlist = range(2001,2015)

    for i in range(len(lossyear)):
        tmp_props = lossyear[i]['properties']

        tmp_total = 0
        if 'loss_none' in tmp_props:
            tmp_total += tmp_props['loss_none']

        for j in range(2001,2015):
            jstr = 'loss_' + str(j)
            if jstr in tmp_props:
                tmp_total += tmp_props[jstr]


        for j in range(2001,2015):
            jstr = 'loss_' + str(j)
            
            tmp_index = int(jstr[-2:]) - 1

            if jstr in tmp_props:
                tmp_percent_lossyear = 100 * tmp_props[jstr] / float(tmp_total)
                percent_lossyear_data[tmp_index].append(roundxy(tmp_percent_lossyear))
            else:
                percent_lossyear_data[tmp_index].append(0)       




    # ltdr_data = [i['properties']['mean'] for i in ltdr]



    tmp_df = pd.DataFrame(tmp_data)
    tmp_df.rename(columns = {'mean':'tc00_e'}, inplace=True)

    tmp_df['tc00_e'] = tmp_df['tc00_e'].apply(roundxy)
    tmp_df['per_loss'] = percent_loss_data
    tmp_df['per_gain'] = percent_gain_data


    # add loss year data
    for i in range(len(percent_lossyear_data)):
        istr = str(i+1)

        if len(istr) == 1:
            year = "200" + istr
        else:
            year = "20" + istr

        tmp_df['per_loss_'+year] = percent_lossyear_data[i]



    # tmp_df['lnyx_2000e'] = ltdr_data

    #tmp_df['ad_extract'].fillna('NA', inplace=True)
    tmp_df.to_csv(base_dir + "/analysis_data/hansen_extracts/"+rid+"_extract.csv", sep=",", encoding="utf-8", index=False)

