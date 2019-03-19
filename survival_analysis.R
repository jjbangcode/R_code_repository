
library(survival)
c_dat = as.matrix(read.table('/home/starjjbang/Storage/part1/breast-ASV/data/070921/ASV/clinical/clinical_105-1.txt',sep='\t'))
c_dat[which(c_dat[,5]=="Dead"),5] = 1
c_dat[which(c_dat[,5]=="Alive"),5] = 0
colnames(c_dat) = c("id","age","stage","dat.S","sur") 
c_dat[which(c_dat[,3]=="Stage IV"),3] = 9
c_dat[which(c_dat[,3]=="Stage IIA"),3] = 8
c_dat[which(c_dat[,3]=="Stage IA"),3] = 7
c_dat[which(c_dat[,3]=="Stage IIB"),3] = 6
c_dat[which(c_dat[,3]=="Stage IIIC"),3] = 5
c_dat[which(c_dat[,3]=="Stage IB"),3] = 4
c_dat[which(c_dat[,3]=="Stage I"),3] = 3
c_dat[which(c_dat[,3]=="Stage IIIA"),3] = 2
c_dat[which(c_dat[,3]=="Stage IIIB"),3] = 1
c_dat[which(c_dat[,3]==NA),3] = 0


#somatic mutation survival
g_dat = as.matrix(read.table('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/6_somatic_)mutation_indi.txt',sep='\t'))
a_lab = unique(g_dat[,1])
pv = c()
for( i in 1:length(a_lab))
{
  tmp = rep(0,length(c_dat[,1]))
  names(tmp) = c_dat[,1]
  s_lab = unique(g_dat[which(g_dat[,1]==a_lab[i]),2])
  tmp[which(names(tmp)==s_lab)]=1
  
  dat2 = c_dat[,2:ncol(c_dat)]
  colnames(dat2) = c('age','stage','dat.S','sur')
  rownames(dat2) = c_dat[,1]
  
  dat3 = as.data.frame(cbind(dat2,tmp))
  colnames(dat3) = c('age','stage','dat.S','sur','exist')
  t = as.numeric(as.character(dat3$dat.S))
  st = as.numeric(as.character(dat3$sur))
  age = as.numeric(as.character(dat3$age))
  sc = as.numeric(as.character(dat3$exist))
  stage = as.character(dat3$stage)
  fit = coxph(Surv(t,st)~ age + sc, method="efron",control = coxph.control(iter.max = 50))
  p.val = summary(fit)$logtest[3]
  pv = rbind(pv,c(a_lab[i],p.val,length(s_lab)))
  print(i)
}

write.table(pv,'/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/20190319_somatic_survival_result.txt'
            ,sep='\t',quote=FALSE,row.names = F)

#SAAV survival
p_dat = as.matrix(read.table('/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/5_SAAV_indi.txt',sep='\t',header = T))
a_lab = unique(p_dat[,1])

p_pv = c()
for(i in 1:length(a_lab))
{
  tmp = rep(0, ncol(p_dat)-1)
  names(tmp) = unlist(lapply(colnames(p_dat)[2:ncol(p_dat)], function(x) gsub('\\.','-',x)))
  lab_tmp = unique(p_dat[which(p_dat[,1] == a_lab[i]),])[2:ncol(p_dat)]
  tmp[which(is.na(lab_tmp)==FALSE)]=1
  tmp = as.matrix(tmp)
  
  dat2 = c_dat[,2:ncol(c_dat)]
  colnames(dat2) = c('age','stage','dat.S','sur')
  rownames(dat2) = c_dat[,1]
  dat3 = merge(dat2,tmp, by='row.names')
  colnames(dat3) = c('id','age','stage','dat.S','sur','exist')
  t = as.numeric(as.character(dat3$dat.S))
  st = as.numeric(as.character(dat3$sur))
  age = as.numeric(as.character(dat3$age))
  sc = as.numeric(as.character(dat3$exist))
  stage = as.character(dat3$stage)
  fit = coxph(Surv(t,st)~ age + sc, method="efron",control = coxph.control(iter.max = 50))
  p.val = summary(fit)$logtest[3]
  p_pv = rbind(p_pv,c(a_lab[i],p.val,length(which(is.na(lab_tmp)==FALSE))))  
  print(i)
}
write.table(pv,'/home/starjjbang/Storage/part1/SAAVpedia/test_result/20190312_SNV_val/20190319_SAAV_survival_result.txt'
            ,sep='\t',quote=FALSE,row.names = F)

