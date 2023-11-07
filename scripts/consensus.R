library(stringr)

files <- list.files(path = './analysis/consensus')

counter <- 0
cons_orig <- list()
for(f in files){
  counter <- counter+1
  print(paste0('creating consensus for ',str_replace(f,'.csv',''),' [',counter,' of ',length(files),']'))
  df <- read.csv(paste0('./analysis/consensus/',f),header = TRUE)
  df <- df[,-1]
  consensus <- c()
  for(pos in 1:dim(df)[1]){
    aa = colnames(df)[which(df[pos,] == max(df[pos,],na.rm = TRUE))]
    if(length(aa) > 1){
      aa = sample(aa,1)
      consensus <- append(consensus , aa)
    }else{
      if(length(aa) == 1){
        if(aa == 'X.'){
          consensus <- append(consensus , '-')
        }else{
          consensus <- append(consensus , aa)
        }
      }
    }
  }
  cons_aligned = substr(paste(consensus,collapse = ''),1069,1184) # region of interest
  co <- str_replace_all(cons_aligned,'-','')
  cons_orig <- append(cons_orig,co)
}

finaltab <- t(data.frame(cons_orig))
rownames(finaltab) <- str_replace_all(files,'.csv','')
colnames(finaltab) <- 'consensus'
write.csv(finaltab,'./analysis/consensus_tab.csv',quote = F)
