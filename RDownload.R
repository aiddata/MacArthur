#Download Script

downlad_data <- function(mDir)
{
  #Flag for if data was downloaded or not
  new_data = FALSE
  
  if(mDir == FALSE)
  {
    this.dir <- dirname(parent.frame(2)$ofile)
    mDir = this.dir
  }
  
  ftp_url <- "ftp://reu:reudata@ftp.aiddata.wm.edu/REU/MacArthur/analysis_data/"
  dir_names = getURL(ftp_url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  dirs_char <- strsplit(dir_names, "\r*\n")[[1]]
  
  dirs_int <- as.numeric(dirs_char)
  dirs_order <- as.character(dirs_int[order(-dirs_int)])
  
  dir_active = dirs_order[1]
  
  active_dir_path = paste(mDir, "/analysis_data/", dir_active, sep='')
  
  # make new dir
  dir.create(active_dir_path, showWarnings = FALSE)
  
  # get files
  active_url <- paste(ftp_url, dir_active, '/', sep='')
  file_str = getURL(active_url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  file_names <- strsplit(file_str, "\r*\n")[[1]]
  
  
  # download
  for (i in 1:length(file_names)) {
    print(i)
    
    
    src_file <- paste(active_url, file_names[[i]], sep="")
    print(src_file)
    
    dst_file <- paste(active_dir_path, '/', file_names[[i]], sep="")
    
    ftype <- length(strsplit(src_file, "[.]")[[1]])-3
    
    if (!file.exists(dst_file) & (ftype > 1)) 
    {
      download.file(src_file, destfile=dst_file, method="libcurl")
    }
    if(!file.exists(dst_file) & (ftype <= 1))
    {
      dir.create(paste(active_dir_path,"/",file_names[[i]], "/", sep=""), showWarnings = FALSE)
      folder_url <- paste(ftp_url, dir_active, '/',file_names[[i]],'/', sep='')
      fold_url = getURL(folder_url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
      fold_file_names <- strsplit(fold_url, "\r*\n")[[1]]
      for(j in 1:length(fold_file_names))
      {
        src_filef <- paste(folder_url, fold_file_names[[j]], sep="")
        dst_filef <- paste(active_dir_path, '/', file_names[[i]], '/', fold_file_names[[j]], sep="")
        download.file(src_filef, destfile=dst_filef, method="libcurl")
      }
    }
  }
  return(active_dir_path)
}