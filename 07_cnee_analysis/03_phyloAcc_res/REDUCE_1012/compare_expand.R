prefix1 = "../NEW_CNEE/reduced/result_2018_1012_gap/subseq_cons"

# combine three runs
options(scipen = 10)
score = matrix(0, 0, 9)
record_run = matrix(1,284001,3)  # each run is the best for each model
for(i in 0:56)
{
  str = i * 5000 + 1 ; ed = (i+1) * 5000
  if(i == 56) ed = 284001
  for(k in 1:3)
  {
    f2 = read.table(paste0(prefix1, "_",str,"_", ed, "_v",k,"_elem_lik.txt"), header=T)
    #print(sum(is.na(f1$logBF1)))
    if(k==1){ f1 = f2 }else{
      for(j in 3:5){
      ind = which(f1[,j] < f2[,j])
      f1[ind,c(j, j+5)]  = f2[ind,c(j, j+5)]
      record_run[str+ind-1, j-2] = k
      }
    }
  }
  f1[,1] = f1[,1] + i * 5000
  score <- rbind(score, f1)
}
score$logBF1 = score$loglik_Acc - score$loglik_Null
score$logBF2 = score$loglik_Acc - score$loglik_Full
write.table(score, "gap/Combined_elem_lik.txt", quote=F, sep="\t", row.names = F)



## combine three runs
for(m in 0:2)
{
  all_Z <- matrix(0, 0, 377) #341, 345,377
  for(i in 0:56)
  {
    str = i * 5000 + 1 ; ed = (i+1) * 5000
    if(i == 56) ed = 284001
    tmplate <- data.frame(matrix(0, ed - str + 1, ncol(all_Z)))
    for(k in 1:3)
    {
      Z <- read.table(paste0(prefix1, "_",str,"_", ed, "_v",k,"_rate_postZ_M",m,".txt"),header=T, row.names = 1) #
      tmplate1 <- data.frame(matrix(0, ed - str + 1, ncol(Z)))
      rownames(tmplate1) <- 0:(ed - str)
      tmplate1[rownames(Z),] <- Z
      ind = which(record_run[str:ed,m+1] == k)
      tmplate[ind,] <- tmplate1[ind,]
    }
    rownames(tmplate) <- (str-1) : (ed-1)
    all_Z <- rbind(all_Z, tmplate)
  }
  
  colnames(all_Z) <- colnames(Z)
  rownames(all_Z) = score$ID
  write.table(all_Z, paste0("gap/Combined_post_Z_M",m,".txt"), quote=F, sep="\t", row.names = T) 
}
