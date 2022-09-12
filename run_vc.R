##############################################################################
# Input the response data and 
# generate the response matrix ygt for vc model fitting
#                 n: number of observations
#		  p: number of predictors
#		  T: number of time points sampled
##############################################################################
rm(list = ls())

rootDir = '/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/3T/3T_2020/SVC/neworder_3T/'# where the data is
fname = paste(rootDir,'/mem_ANTN_bin10_norm.csv', sep = '') # create full path for input data
setwd(paste(rootDir,sep = '//')) # set it at where the vc scripts are

response=read.csv(fname,header=TRUE)#1st column:memory; 2nd Column: bin

qq=response[1]
response=subset(response,!is.na(qq))

dim(response); # check data dimensions 
n= dim(response)[1]; # number of observations (~ nSubj)
tVec=unique(response$bin) #remove repeated 
tVec=sort(tVec)
T=length(tVec); 

for (i in 1:dim(response)[1]){
  response[i,'t']=which (tVec==response[i,'bin']) # 
}
colnames(response)=c("memory", "bin","t")

###### Obtain n*T matrix ygt and set the missing value of Y as 0
ygt = array(0,dim=c(n,T))


for (nI in 1:n){
  ygt[nI,response[nI,3]]=response[nI,1] #
  # ygt[nI,response[nI,'t']]=response[nI,'amytau']
}


image(ygt,col= grey(seq(0,1,length = 256)), axes = F, xlab = 'SID',ylab = 'amytau')
axis(2, at = seq(0,1,length.out  = T), labels = tVec)
axis(1, at = seq(0,1,length.out = n), labels = 1:n)


##############################################################################
# Input the predictor variables and 
# generate the predictor matrix X for vc model fitting
#                 n: number of observations
#		  p: number of predictors
#		  T: number of time points sampled
#############################################################################
fname_predic = paste(rootDir,'/brainscores_z_gender_edu_APOE_ICV_ANTN.csv',sep = '')
predic<-read.csv(fname_predic,header=TRUE)#

predic=subset(predic,!is.na(qq))

dim(predic)
p= dim(predic)[2] 

#Standardize predic matrix, define custom function
stand_matrix<- function (mat){
  
  normX = rep(0,dim(mat)[2])
  meanX = rep(0,dim(mat)[2])
  
  for (pI in 1:dim(mat)[2]) {
    
    
    meanX[pI] = mean(mat[,pI]);
    normX[pI] = sqrt( sum( (mat[,pI]-meanX[pI])^2 ) )
    if (normX[pI]==0) { normX[pI]= 1 }
    mat[,pI] = (mat[,pI] - meanX[pI])/normX[pI]
  }
  return(mat)
}

predic_std<-stand_matrix(predic)


# Create standardized predictor matrix X based on predic_std (for missing data we set x as 0)
X_std = array(0, dim=c(n,p,T));
dim(X_std)

for (nI in 1:n){
  for (pI in 1:p){
    #X_std[nI,pI,response$t[nI]]<- predic_std[nI,pI]
    X_std[nI,pI,response$t[nI]]<- predic[nI,pI]
  }
}


#############################################################
# run the main function to estimate beta(t) by running vc_sparse, default NO INTERCEPT
########################################################
W = matrix(0,p,p) # weighted matrix for predictors , but not used in this project

source(paste("/mnt/isilon/CSC4/HelenZhouLab/HZLHD2/Data7/Members/Qianxing/ADNI_NgKP/3T/SVC/vc.R",sep = '')) # 

################################################################

#####using for loop to repeat the VC model fitting, default NO INTERCEPT (need to change in aging.R vc_sparse)
sim_no=100
beta_result=array(0,dim=c(p,T,sim_no))
for (i in 1:sim_no){
  set.seed(i+1000)
  res_std=vc_sparse(simLen=1,n=n, T=T, p=p,  X=X_std,ygt=ygt,W=W,tVec=tVec,
                    do.alpha2=FALSE, do.alpha3=TRUE)
  beta_result[,,i]=res_std$beta
  
  # Display for loop count
  cat(i,'  ')
  if (i==50){ cat('\n') }
  if (i==100){ cat('\n') }
}

save(list=ls(), file=paste(rootDir, "mem_brainscores_gender_edu_APOE_ICV_ANTN_bin10_norm", Sys.Date() , ".RData", sep = ''))



beta_sum<-rep(0,length=sim_no) # count how many predictors are kept in each vc loop
for (i in 1:sim_no){
  beta_sum[i]<- length(which(rowSums(beta_result[,,i])!=0))
  #cat(which(rowSums(beta_result[,,i])!=0))
  #cat("\n")}
}

#######

############Dispaly all the variables which are selected
# display index
for (i in 1:sim_no){
  cat(which(rowSums(beta_result[,,i])!=0))
  cat("\n")}

# display variable name
#for (i in 1:sim_no){
#  cat('sim_no = ',i,'\n',sep = '')
#  print(colnames(predic)[which(rowSums(beta_result[,,i])!=0)])
#  cat("\n")}

## Track the number of times a predictor is selected
vCount = array(0,dim = c(p,sim_no))
for (i in 1:sim_no){
  vCount[which(rowSums(beta_result[,,i])!=0),i] = 1
}

which(rowSums(vCount) > 0)





selection_criterion = sim_no # how many times a predictor were selected across all vc iterations. max = sim_no
var_idx = which(rowSums(vCount) == selection_criterion) #==
chart_title = colnames(predic)[var_idx]
var_idx

selection_criterion = 100#sim_no # how many times a predictor were selected across all vc iterations. max = sim_no
var_idx = which(rowSums(vCount) >= selection_criterion) #==
chart_title = colnames(predic)[var_idx]
var_idx

#compute the mean and SD of beta_result
beta_result_mean=array(0,dim=c(p,T))
beta_result_SD=array(0,dim=c(p,T))


for (i in 1: p){
  for (j in 1: T){
    
    beta_result_mean[i,j]= mean(beta_result[i,j,])
    beta_result_SD[i,j]  = sd(beta_result[i,j,])
  }
}


length(var_idx)
par(mfrow=c(3,4))

for (vI in 1: length(var_idx)){ 
  chart_title= colnames(predic)[var_idx[vI]]
  
  plot(tVec, beta_result_mean[var_idx[vI],],xlab="bin",ylab="Beta",type='l')#)
  lines(tVec, beta_result_mean[var_idx[vI],]+2*beta_result_SD[var_idx[vI],]/sqrt(sim_no), col="red",lty=3)
  lines(tVec, beta_result_mean[var_idx[vI],]-2*beta_result_SD[var_idx[vI],]/sqrt(sim_no), col="red",lty=3)
  
  lines(tVec, array(0,dim=c(1,length(tVec))),add=TRUE)
  # lines(array(33,dim=c(1,round(max(beta_result_mean[var_idx[vI],]))-round(min(beta_result_mean[var_idx[vI],]))+1)),round(min(beta_result_mean[var_idx[vI],])):round(max(beta_result_mean[var_idx[vI],])),add=TRUE)
  # lines(array(69,dim=c(1,round(max(beta_result_mean[var_idx[vI],]))-round(min(beta_result_mean[var_idx[vI],]))+1)),round(min(beta_result_mean[var_idx[vI],])):round(max(beta_result_mean[var_idx[vI],])),add=TRUE)
  
  
  title(chart_title,line=3)
}

#write.csv(beta_result_mean, file = "mem_brainscores_gender_edu_APOE_ICV_ANTN_bin10_norm_ridge_beta_mean.csv", append = FALSE, quote = TRUE, sep = " ",
 #         eol = "\n", na = "NA", dec = ".", row.names = TRUE,
 #         col.names = TRUE, qmethod = c("escape", "double"),
 #         fileEncoding = "")
#write.csv(beta_result_SD, file = "mem_brainscores_gender_edu_APOE_ICV_ANTN_bin10_norm_ridge_beta_SD.csv", append = FALSE, quote = TRUE, sep = " ",
  #        eol = "\n", na = "NA", dec = ".", row.names = TRUE,
  #        col.names = TRUE, qmethod = c("escape", "double"),
  #        fileEncoding = "")


#dev.copy2eps(file="2variables_SD_plot.eps")
dev.off()












