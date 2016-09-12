#!/bin/bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

# https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.2.html
# https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/treecover2000.txt

# region list containing region ids
rlist=(sea sa africa)
# rlist=(sea)

for rid in ${rlist[@]}; do

    echo ${rid}

    extent_grep=$(ogrinfo ../shps/${rid}/${rid}_grid.shp ${rid}_grid | grep Extent)

    if [[ ${extent_grep} != Extent:* ]]; then
        echo 'No extent data found for '${rid}'_grid.shp'
        exit 1
    fi

    mkdir -p downloads/{treecover2000,loss,gain,lossyear}/${rid}_tiles

    echo $extent_grep

    # examples of grep output, format and needed values
    # Extent: (92.189278, 5.593261) - (109.489278, 28.543261)
    # Extent: (LL lon, LL lat) - (UR lon, UR lat)
    # Extent: (90, 0) - (110, 30)

    extent_string=`echo $extent_grep | sed s/'Extent: '//g | sed s:' - ':' ':g | sed s:[\),\(]::g`
    echo $extent_string

    cmd='python bbox_to_filenames.py '${extent_string}
    # echo $cmd

    url_string=$(eval $cmd)
    # echo $url_string

    url_array=($url_string)

    for i in ${url_array[@]}; do
        echo $i

        # wget tiles

        wget -c -N -P downloads/treecover2000/${rid}_tiles ${i}

        j=`echo $i | sed s/'treecover2000'/'loss'/g`
        wget -c -N -P downloads/loss/${rid}_tiles ${j}

        k=`echo $i | sed s/'treecover2000'/'gain'/g`
        wget -c -N -P downloads/gain/${rid}_tiles ${k}

        m=`echo $i | sed s/'treecover2000'/'lossyear'/g`
        wget -c -N -P downloads/lossyear/${rid}_tiles ${m}

    done

    mkdir -p data/{treecover2000,loss,gain,lossyear}


    # mosaic tiles

    gdal_merge.py -of GTiff -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=YES downloads/treecover2000/${rid}_tiles/*.tif -o data/treecover2000/${rid}_mosaic.tif
    gdal_merge.py -of GTiff -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=YES downloads/loss/${rid}_tiles/*.tif -o data/loss/${rid}_mosaic.tif
    gdal_merge.py -of GTiff -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=YES downloads/gain/${rid}_tiles/*.tif -o data/gain/${rid}_mosaic.tif
    gdal_merge.py -of GTiff -co COMPRESS=LZW -co TILED=YES -co BIGTIFF=YES downloads/lossyear/${rid}_tiles/*.tif -o data/lossyear/${rid}_mosaic.tif


done

