rm(list=ls())

library(rstudioapi)
current_path <-rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

library(kknn)
library(SNFtool)
library(MASS)
library(flexclust)
library(coda)
library(ISLR)
library(pROC)
library(ROCR)
library(ggplot2)
library(dplyr)
library(reshape)
library(MCMCpack)

#---------------------------------------------------------------------------------------------------------------
# Basic Settings // Data Loading
#---------------------------------------------------------------------------------------------------------------
nitem = 72;
nschool = 62;
ndyad = choose(nitem, 2);
ndim = 2;
ncat = 1;
nmax = 81;
niter = 2500;
count = scan("DATA/count.txt")
option = 1

# Load Dataset and Make Adjacency Matrices
itemsum = matrix(NA,nschool,nitem)
sampsum = matrix(NA,nschool,nmax)
Y = array(0, dim=c(nschool,nitem,nmax,nmax))
U = array(0, dim=c(nschool,nmax,nitem,nitem))
for(a in 1:nschool){
  if(a < 10){
    vname = paste("item0",a,sep="")
    fname = paste("DATA/item0",a,".txt",sep="")
  }
  else{
    vname = paste("item",a,sep="")
    fname = paste("DATA/item",a,".txt",sep="")
  }
  temp = matrix(scan(fname),ncol=nitem,byrow=TRUE)
  assign(vname,temp)
  
  for(k in 1:nitem){
    for(i in 2:count[a]){
      for(j in 1:(i-1)){
        Y[a,k,i,j] = Y[a,k,j,i] = temp[i,k] * temp[j,k]
      }
    }
  }
  for(k in 1:count[a]){
    for(i in 2:nitem){
      for(j in 1:(i-1)){
        U[a,k,i,j] = U[a,k,j,i] = temp[k,i] * temp[k,j]
      }
    }
  }
  itemsum[a,] = colSums(temp)
  sampsum[a,(1:count[a])] = rowSums(temp)
}

itemprob = matrix(0,nschool,nitem)
for(a in 1:nschool) itemprob[a,] = itemsum[a,]/count[a]
itemjump = ceiling(itemprob * 10)
write.table(itemjump, "DATA/jumpitem.txt", row.names=FALSE, col.names=FALSE)

sim_beta  = array(NA,dim=c(nschool,niter,nitem))
sim_theta = array(NA,dim=c(nschool,niter,nmax))
sim_sigz  = matrix(NA,niter,nschool)
sim_z = array(NA,dim=c(nschool,niter,nmax*ndim))
sim_i = array(NA,dim=c(nschool,niter,nitem*ndim))
theta_result = rep(NA,sum(count))
a = 1:niter; index = 0;
for(i in 1:nschool){
  if(i < 10){
    fname1 = paste("RESULT/sim0",i,"_b1.log",sep="")
    fname2 = paste("RESULT/sim0",i,"_t1.log",sep="")
    fname3 = paste("RESULT/sim0",i,"_h1.log",sep="")
    fname4 = paste("RESULT/sim0",i,"_z1.log",sep="")
    fname5 = paste("RESULT/sim0",i,"_i1.log",sep="")
  }
  else{
    fname1 = paste("RESULT/sim",i,"_b1.log",sep="")
    fname2 = paste("RESULT/sim",i,"_t1.log",sep="")
    fname3 = paste("RESULT/sim",i,"_h1.log",sep="")
    fname4 = paste("RESULT/sim",i,"_z1.log",sep="")
    fname5 = paste("RESULT/sim",i,"_i1.log",sep="")
  }
  temp = matrix(scan(fname1),ncol=nitem,byrow=TRUE)
  sim_beta[i,,] = temp
  
  temp = matrix(scan(fname2),ncol=count[i],byrow=TRUE)
  sim_theta[i,,1:count[i]] = temp
  theta_result[(index+1):(index+count[i])] = apply(temp,2,mean)
  index = index + count[i]
  
  temp = matrix(scan(fname4),ncol=count[i]*ndim,byrow=TRUE)
  sim_z_dist = array(NA,dim=c(niter,count[i],count[i]))
  sim_z[i,,1:(count[i]*ndim)] = temp
  for(iter in 1:niter){
    temp_mat = matrix(temp[iter,],count[i],ndim,byrow=TRUE)
    sim_z_dist[iter,,] = as.matrix(dist(temp_mat))
  }
  
  temp = matrix(scan(fname5),ncol=nitem*ndim,byrow=TRUE)
  sim_i_dist = array(NA,dim=c(niter,nitem,nitem))
  sim_i[i,,1:(nitem*ndim)] = temp
  for(iter in 1:niter){
    temp_mat = matrix(temp[iter,],nitem,ndim,byrow=TRUE)
    sim_i_dist[iter,,] = as.matrix(dist(temp_mat))
  }
  
  # Get $\sigma_{\epsilon_k}^2$ for each school
  temp = scan(fname3)
  sim_sigz[,i] = temp
}

item_index = matrix(NA,choose(nitem,2),2)
temp = 0
for(i in 2:nitem){
  for(j in 1:(i-1)){
    temp = temp + 1
    item_index[temp,1] = i
    item_index[temp,2] = j
  }
}

school_index = matrix(NA,choose(nschool,2),2)
temp = 0
for(i in 2:nschool){
  for(j in 1:(i-1)){
    temp = temp + 1
    school_index[temp,1] = i
    school_index[temp,2] = j
  }
}

sim_gamma  = matrix(scan("RESULT/sim_g1.log"),ncol=nitem,byrow=TRUE)
sim_varphi = matrix(scan("RESULT/sim_p1.log"),ncol=nitem,byrow=TRUE)
sim_delta  = matrix(scan("RESULT/sim_l1.log"),ncol=ndyad*ncat,byrow=TRUE)
sim_tau    = matrix(scan("RESULT/sim_u1.log"),ncol=ndyad*ncat,byrow=TRUE)

#---------------------------------------------------------------------------------------------------------------
# Mean and 95% HPD intervals for $\gamma$
# Figure 2(b) of the manuscript
#---------------------------------------------------------------------------------------------------------------
gamma_mean = apply(sim_gamma,2,mean)
gamma_mcmc = mcmc(sim_gamma)
gamma_result = matrix(NA,length(gamma_mean),3)
gamma_result[,1] = gamma_mean
gamma_result[,2:3] = HPDinterval(gamma_mcmc,prob=0.95)
round(gamma_result,4)
min(gamma_result)
max(gamma_result)

pdf("PLOT/hpd_gamma.pdf")
x.axis = c(1:nitem)
y.axis = x.axis
plot(x.axis,y.axis,ylim=c(-3.75,2.50),xlab="item",ylab=expression(gamma),type="n")
for(i in 1:nitem){
  lines(c(i,i),gamma_result[i,2:3],lwd=2)
  points(i,gamma_result[i,1],pch=20,col=4)
} 
for(i in 1:7) abline(v=i*10, lty=3, col=2)
dev.off()


temp = matrix(scan("RESULT/sim_g1.log"),ncol = 72,byrow=TRUE)
temp_dist <- matrix(0, nrow = 2500 * 72, ncol = 1)
for(j in 1:72){
  for(i in 1:2500){
    temp_dist[i + (j-1) * 2500] <- temp[2*i, j] - temp[2*i - 1, j]
  }
}

pdf("GGPLOT/hpd_gamma_diff.pdf")
gamma_result2 <- data.frame(value = temp_dist, item = rep(1:72, each = 2500))
gamma_result3 <- data.frame(gamma_result2 %>% group_by(item) %>% summarise(min = min(value), max = max(value), mean = mean(value)))
gg <- ggplot(gamma_result3)+
  geom_linerange(aes(x=item, ymin=min, ymax=max), size = 1) +
  geom_point(aes(x=item, y=mean), color = "blue") +
  geom_vline(xintercept=round(seq(from = 0, to = 72, by = 10)), linetype = 'dotted', color='red', size = 0.7) +
  scale_x_continuous(breaks = round(seq(from = 0, to = 72, by = 10))) +
  xlab("Item") + ylab(expression(gamma))
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Visualization of the item latent positions using from estimated item distance matrix at the upper level.
# 1. Heatmap for estimated $\mu$
# 2. Reconstruction of the item latent positions using from $\mu$
#---------------------------------------------------------------------------------------------------------------
temp = matrix(scan("RESULT/sum_l1.log"),ncol=ndyad,byrow=TRUE)
temp1 = temp[1,1:ndyad];
temp2 = temp[3,1:ndyad];

# Traditional School System
idist_control = matrix(0,nitem,nitem);
logidist_control = matrix(0,nitem,nitem);
index = 0;
for(i in 2:nitem){
  for(j in 1:(i-1)){
    index = index + 1
    idist_control[i,j] = exp(temp1[index])
    logidist_control[i,j] = temp1[index]
  }
}
idist_control = idist_control + t(idist_control)
logidist_control = logidist_control + t(logidist_control)

if(option == 1){
  ntrad = isoMDS(idist_control,k=2)
  ntrad_x = ntrad$points[,1]; ntrad_y = ntrad$points[,2];
  write.table(ntrad$points,"RSUM/position_control.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/position_control.txt"),ncol=2,byrow=TRUE)
  ntrad_x = temp[,1]; ntrad_y = temp[,2];
}

# Renovation School System
idist_innovation = matrix(0,nitem,nitem);
logidist_innovation = matrix(0,nitem,nitem);
index = 0;
for(i in 2:nitem){
  for(j in 1:(i-1)){
    index = index + 1
    idist_innovation[i,j] = exp(temp2[index])
    logidist_innovation[i,j] = temp2[index]
  }
}
idist_innovation = idist_innovation + t(idist_innovation)
logidist_innovation = logidist_innovation + t(logidist_innovation)

if(option == 1){
  ninno = isoMDS(idist_innovation,k=2)
  #ntrad = cmdscale(idist,eig=TRUE,k=2)
  ninno_x = ninno$points[,1]; ninno_y = ninno$points[,2];
  write.table(ninno$points,"RSUM/position_innovation.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/position_innovation.txt"),ncol=2,byrow=TRUE)
  ninno_x = temp[,1]; ninno_y = temp[,2];
}

#rotation
#---------------------------------------------------------------------------------------------------------------
# Figure 5 in the Manuscript
# Reconstruction of the item latent positions using from $\mu$
#---------------------------------------------------------------------------------------------------------------
temp1 = ntrad$points
temp2 = ninno$points

pdf("GGPLOT/item_ls_innovation.pdf")
temp_ls_inno <- data.frame(temp2)
gg <- ggplot() + 
  geom_text(data = temp_ls_inno, aes(x = X1, y = X2, label=1:72), size = 4.5, fontface = "bold") +
  xlab(NULL) + ylab(NULL)
print(gg)
dev.off()

pdf("GGPLOT/item_ls_control.pdf")
temp_ls_cont <- data.frame(temp1)
gg <- ggplot() + 
  geom_text(data = temp_ls_cont, aes(x = X1, y = X2, label=1:72), size = 4.5, fontface = "bold") +
  xlab(NULL) + ylab(NULL)
print(gg)
dev.off()

rot_position <- procrustes(as.matrix(temp_ls_inno), as.matrix(temp_ls_cont))$X.new
temp_ls_inno2 = data.frame(rot_position)

pdf("GGPLOT/item_ls_innovation_rot.pdf")
gg <- ggplot() + 
  geom_text(data = temp_ls_inno2, aes(x = X1, y = X2, label=1:72), size = 4.5, fontface = "bold") +
  xlab(NULL) + ylab(NULL)
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Figure 6 of the Manuscript
# Heatmap for estimated $\mu$
#---------------------------------------------------------------------------------------------------------------
#control heatmap
hist_temp <- as.matrix(exp(-dist(temp_ls_cont)))
dim(hist_temp)
pdf("GGPLOT/heatmap_control.pdf");
heatmap(hist_temp, symm=TRUE, keep.dendro=FALSE, cexRow=0.5, cexCol=0.5)
dev.off();

# innovation heatmap
hist_temp <- as.matrix(exp(-dist(temp_ls_inno)))
dim(hist_temp)
pdf("GGPLOT/heatmap_innovation.pdf");
heatmap(hist_temp, symm=TRUE, keep.dendro=FALSE, cexRow=0.5, cexCol=0.5)
dev.off();

#---------------------------------------------------------------------------------------------------------------
# Save Computation Result
#---------------------------------------------------------------------------------------------------------------
save.image("GEPS_REVERSE_MODEL2.RData")
