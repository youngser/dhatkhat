## access Google Drive
library(googledrive)
folder_url <- "https://drive.google.com/open?id=1zBuIMdpNQFoWZvh7JOepswwFCq2l9Wco"
folder <- drive_get(as_id(folder_url))
files <- drive_ls(folder)
f1 <- drive_download(files[1,])