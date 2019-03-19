#somatic mutation
g_dat = as.matrix(read.csv('/home/starjjbang/Storage/part1/SAAVpedia/data/breast_cancer_genome/somatic.csv',sep=','))

g_dat1 =c()
for(i in 1:nrow(g_dat))
{
  if( g_dat[i,6]=='SNP')
  {
    k = ''
    if(g_dat[i,7] == g_dat[i,8])
    {
      k = paste(unlist(lapply(g_dat[i,c(1,2,7,9)], function(x) trimws(x))),collapse = ';' )
    }
    else
    {
      k = paste(unlist(lapply(g_dat[i,c(1,2,7,8)], function(x) trimws(x))),collapse = ';' )
    }
    g_dat1 = rbind(g_dat1, c(k,paste(unlist(strsplit(g_dat[i,12],'-'))[2:3],collapse = '-')))
  }
}
write.table(g_dat1,'/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/6_somatic_)mutation_indi.txt',sep='\t',quote=FALSE, row.names = FALSE, col.names = F)

#SAAV
p_dat = as.matrix(read.csv('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/2_SAAV_pos.txt',sep='\t',header=F))
colnames(p_dat) = p_dat[1,]
p_dat = p_dat[2:nrow(p_dat),]
p_tmp = c()
for(i in 1:nrow(p_dat))
{
  k_make = lapply(unlist(strsplit(p_dat[i,4],',')), function(x)  paste(p_dat[i,1],p_dat[i,2],p_dat[i,3],x,sep = ';'))
  for( j in 1:length(k_make))
  {
    p_tmp = rbind(p_tmp,c(k_make[j],p_dat[i,13:ncol(p_dat)] ))
  }
}
head(p_tmp)
write.table(p_tmp,'/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/5_SAAV_indi.txt',sep='\t',quote=FALSE, row.names = FALSE, col.names = T)

#paring individisual 
g_dat = as.matrix(read.table('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/6_somatic_)mutation_indi.txt',sep='\t'))
p_dat = as.matrix(read.table('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/5_SAAV_indi.txt',sep='\t',header = T))

map_dat <- function(lab_dat, input_dat)
{
  tmp = rep(0,length(lab_dat))
  tmp[which(lab_dat %in% input_dat)]=1
  return(tmp)
}

for(i in 2:ncol(p_dat))
{
  p_tmp = p_dat[,c(1,i)]
  p_tmp1 = p_tmp[which(is.na(p_tmp[,2])==FALSE),]
  g_tmp = g_dat[which(g_dat[,2]==gsub('\\.','-',colnames(p_dat)[i])),]
  lab =unique(c(g_tmp[,1], t_tmp[,1]))
  a = cbind(lab,map_dat(lab,g_tmp[,1]),map_dat(lab,p_tmp[,1]))
  fn = paste('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190319_indi_SAAV_SNV/',gsub('\\.','-',colnames(p_dat)[i]),'.txt',sep='')
  write.table(a,fn,quote=FALSE,row.names = FALSE,col.names = FALSE,sep='\t')
  print(i)
}

