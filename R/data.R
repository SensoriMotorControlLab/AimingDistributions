
getData <- function() {
  
  dir.create('data/summaries')
  
  Reach::downloadOSFdata(repository='6g3h7',
                         filelist = list('Learners'=c('SummaryFiles.zip')),
                         folder='data/summaries/',
                         overwrite = TRUE,
                         unzip = TRUE,
                         removezips = TRUE,
                         wait=0)
  
}