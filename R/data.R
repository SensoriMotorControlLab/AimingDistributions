
getData <- function() {
  
  if (!dir.exists('data/summaries')) {
    dir.create('data/summaries')
  }
  
  Reach::downloadOSFdata(repository='6g3h7',
                         filelist = list('Learners'=c('SummaryFiles.zip')),
                         folder='data/summaries/',
                         overwrite = TRUE,
                         unzip = TRUE,
                         removezips = TRUE,
                         wait=0)
  
  csvfiles <- list.files('data/summaries/SummaryFiles/', pattern = '.csv')
  for (csvfile in csvfiles) {
    file.rename( from = sprintf('data/summaries/SummaryFiles/%s', csvfile),
                 to = sprintf('data/summaries/%s', csvfile))
  }
  
  # clean up the data directory
  file.remove('data/summaries/SummaryFiles')
  propfiles <- list.files('data/summaries/__MACOSX/SummaryFiles/', all.files=TRUE) # including hidden files
  for (propfile in propfiles) {
    if (propfile %in%  c('.', '..')) { # excluding system directory shortcuts that are not actually files
      next
    } else {
      file.remove(sprintf('data/summaries/__MACOSX/SummaryFiles/%s', propfile))
    }
  }
  file.remove('data/summaries/__MACOSX/SummaryFiles/')
  file.remove('data/summaries/__MACOSX/')
  
}