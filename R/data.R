
getData <- function() {
  
  # if (!dir.exists('data/summaries')) {
  #   dir.create('data/summaries')
  # }
  
  Reach::downloadOSFdata(repository='6g3h7',
                         filelist = list('data'=c('Summaries.zip', 'demographics.csv')),
                         folder='data/',
                         overwrite = TRUE,
                         unzip = TRUE,
                         removezips = FALSE,
                         wait=0)

  
}