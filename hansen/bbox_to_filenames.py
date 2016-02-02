

# cells referenced by top and left
# absolute numbers with NWSE 
# 
# example url
#   https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_treecover2000_10N_100W.tif


import sys
import math
import itertools


ll_lon = float(sys.argv[1])
ll_lat = float(sys.argv[2])
ur_lon = float(sys.argv[3])
ur_lat = float(sys.argv[4])

ll_lon_rnd = int(math.floor(ll_lon / 10) * 10)
ll_lat_rnd = int(math.floor(ll_lat / 10) * 10)
ur_lon_rnd = int(math.ceil(ur_lon / 10) * 10)
ur_lat_rnd = int(math.ceil(ur_lat / 10) * 10)

tops = range(ll_lat_rnd, ur_lat_rnd + 1, 10)[1:]
lefts = range(ll_lon_rnd, ur_lon_rnd, 10)

tops_str = []
lefts_str = []

for i in tops:
    if i < 0:
        out = str(abs(i)) + 'S'
    else:
        out = str(abs(i)) + 'N'

    if len(out) != 3:
        out = '0' + out

    tops_str.append(out)


for i in lefts:
    if i > 0:
        out = str(abs(i)) + 'E'
    else:
        out = str(abs(i)) + 'W'

    if len(out) != 4:
        out = '0' + out

    lefts_str.append(out)


combos = list(itertools.product(tops_str, lefts_str))

# print combos

base = 'https://storage.googleapis.com/earthenginepartners-hansen/GFC2015/Hansen_GFC2015_treecover2000_' #+ '00N' + '_' + '050W' + '.tif'

urls = [base + item[0] + '_' + item[1] + '.tif' for item in combos]

output = ' '.join(urls)

print output

