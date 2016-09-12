
import sys
import os
import time
import pandas as pd

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))



def pandas_fields_monthly_to_yearly(df, field_tuples_list):
    
    if not isinstance(field_tuples_list, list):
        sys.exit("pandas_fields_monthly_to_yearly: invalid field_tuples_list given (must be list)")

    for i in field_tuples_list:
        if not isinstance(i, tuple) and len(i) == 2:
            sys.exit("pandas_fields_monthly_to_yearly: invalid field_tuples_list given (list must contain tuples with 2 items)")

    tmp_df = df.copy(deep=True)

    all_cols = tmp_df.columns.values.tolist()
    
    for i in field_tuples_list:
        
        if i[0].endswith("_"):
            i[0] = i[0][:-1]

        field_cols = [j for j in all_cols if j.startswith(i[0]) and j.endswith(i[1])]

        if len(field_cols) == 0:
            sys.exit("pandas_fields_monthly_to_yearly: no fields found matching starting with "+i[0]+" and ending with "+i[1]+".")


        years = list(set([j.split("_")[1][:4] for j in field_cols]))
        years.sort()

        for y in years:

            # get month fields for year
            year_cols = [j for j in field_cols if j.startswith(i[0]+"_"+y)]
            # print year_cols

            # min
            tmp_df[i[0]+"_"+str(y)+"m"] = tmp_df[year_cols].min(axis=1)

            # mean
            tmp_df[i[0]+"_"+str(y)+"e"] = tmp_df[year_cols].mean(axis=1)
            
            # print tmp_df[i[0]+"_"+str(y)+"e"]
            # max
            tmp_df[i[0]+"_"+str(y)+"x"] = tmp_df[year_cols].max(axis=1)
            
            # drop month fields
            tmp_df.drop(year_cols, axis=1, inplace=True)


    tmp_df.sort_index(axis=1)
    return tmp_df



rlist = ["sea", "sa", "africa"]
# rlist = ["sea"]
# rlist = ["sa"]
# rlist = ["africa"]

for rid in rlist:
    
    print rid


    grid_info_path = base_dir + "/sciclone_data/grid_info/" + rid + "_data.csv"

    hansen_extract_path = base_dir + "/sciclone_data/hansen_extracts/" + rid + "_extract.csv"

    bulk_extract_01_path = base_dir + "/sciclone_data/bulk_extracts/run_01//merge_" + rid + "_grid_edit.csv"

    bulk_extract_02_path = base_dir + "/sciclone_data/bulk_extracts/run_02//merge_" + rid + "_grid.csv"


    grid_info = pd.read_csv(grid_info_path)

    hansen_extract = pd.read_csv(hansen_extract_path)
    # hansen_extract.drop('lnyx_2000e', axis=1, inplace=True)

    bulk_extract_01 = pd.read_csv(bulk_extract_01_path)
    # bulk_cols = bulk_extract_01.columns.values.tolist()
    # cut_cols = [j for j in bulk_cols if j.startswith(('at41', 'pc41'))]
    # cut_bulk_extract = bulk_extract.drop(cut_cols, axis=1)
    # cut_bulk_extract_path = base_dir + "/sciclone_data/bulk_extracts/" + [d for d in os.listdir(base_dir+"/sciclone_data/bulk_extracts") if d.endswith(rid)][0] +"/merge_" + rid + "_grid_edit.csv"
    # cut_bulk_extract.to_csv(cut_bulk_extract_path, index=False, encoding="utf-8")

    bulk_extract_02 = pd.read_csv(bulk_extract_02_path)
    # aggregate year-month fields to year
    agg_bulk_extract_02 = pandas_fields_monthly_to_yearly(bulk_extract_02, [('at41','e'), ('pc41', 'e')])
    # agg_bulk_extract.to_csv(base_dir + "/sciclone_data/merged/test_merge.csv", index=False, encoding="utf-8")


    # merge
    
    tmp_merge = grid_info.merge(hansen_extract, on="ID")

    bulk_merge = bulk_extract_01.merge(agg_bulk_extract_02, on="ID")
    
    final_merge = tmp_merge.merge(bulk_merge, on="ID")

    output_path = base_dir + "/sciclone_data/merged/" + rid + "_merge.csv"
    final_merge.to_csv(output_path, index=False, encoding="utf-8")




# # UNTESTED

# # upload to ftp
# ftp_server = "ftp.aiddata.wm.edu"
# ftp_path = ftp_server + "/REU/MacArthur/analysis_data/" + str(int(time.time()))


# import ftplib

# ftp = ftplib.FTP(ftp_server,'reu','reudata')


# # https://stackoverflow.com/questions/10644608/create-missing-directories-in-ftplib-storbinary

# # Change directories - create if it doesn't exist
# def chdir(dir): 
#     if directory_exists(dir) is False: # (or negate, whatever you prefer for readability)
#         ftp.mkd(dir)
#     ftp.cwd(dir)

# # Check if directory exists (in current location)
# def directory_exists(dir):
#     filelist = []
#     ftp.retrlines('LIST',filelist.append)
#     for f in filelist:
#         if f.split()[-1] == dir and f.upper().startswith('D'):
#             return True
#     return False


# chdir(ftp_path)

# # upload files
# file = open('filename.ext','rb')                  # file to send
# ftp.storbinary('STOR filename.ext', file)     # send the file
# file.close()                                    # close file and FTP
# ftp.quit()



