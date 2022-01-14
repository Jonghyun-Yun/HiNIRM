#---------------------------------------------------------------------------------------------------------------
# Basic Settings // Data Loading
#---------------------------------------------------------------------------------------------------------------
#setwd("~/Documents/BACKUP/MODEL2/RUN1")
rm(list=ls())

library(rstudioapi)
current_path <-rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))

library(kknn);
library(SNFtool);
library(MASS);
library(flexclust);

#---------------------------------------------------------------------------------------------------------------
# Basic Settings // Data Loading
#---------------------------------------------------------------------------------------------------------------
nitem = 72;
nschool = 62;
ndyad = choose(nitem, 2);
ndim = 2;
ncat = 2;
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

#---------------------------------------------------------------------------------------------------------------
# 1. Trace Plots for beta at the lower-level
# 2. Trace Plots for theta at the lower-level
# 3. Trace Plots for the distance between pairs of respondent latent positions at the lower-level
# 4. Trace Plots for the distance between pairs of item latent positions at the lower-level
#---------------------------------------------------------------------------------------------------------------
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
    oname1 = paste("TRACE/beta0",i,".pdf",sep="")
    oname2 = paste("TRACE/theta0",i,".pdf",sep="")
    oname3 = paste("TRACE/resp0",i,".pdf",sep="")
    oname4 = paste("TRACE/item0",i,".pdf",sep="")
  }
  else{
    fname1 = paste("RESULT/sim",i,"_b1.log",sep="")
    fname2 = paste("RESULT/sim",i,"_t1.log",sep="")
    fname3 = paste("RESULT/sim",i,"_h1.log",sep="")
    fname4 = paste("RESULT/sim",i,"_z1.log",sep="")
    fname5 = paste("RESULT/sim",i,"_i1.log",sep="")
    oname1 = paste("TRACE/beta",i,".pdf",sep="")
    oname2 = paste("TRACE/theta",i,".pdf",sep="")
    oname3 = paste("TRACE/resp",i,".pdf",sep="")
    oname4 = paste("TRACE/item",i,".pdf",sep="")
  }
  # Draw the trace plots for $\beta$ for each school
  temp = matrix(scan(fname1),ncol=nitem,byrow=TRUE)
  if(option == 1) {
    pdf(oname1)
    for(j in 1:nitem) plot(a,temp[,j],type="l",main=j)
    dev.off()
  }
  sim_beta[i,,] = temp
  
  # Draw the trace plots for $\theta$ for each school 
  temp = matrix(scan(fname2),ncol=count[i],byrow=TRUE)
  if(option == 1) {pdf(oname2); for(j in 1:count[i]) plot(a,temp[,j],type="l",main=j); dev.off();}
  sim_theta[i,,1:count[i]] = temp
  theta_result[(index+1):(index+count[i])] = apply(temp,2,mean)
  index = index + count[i]
  
  # Draw the trace plots for pairwise distance of respondent latent positions for each school
  temp = matrix(scan(fname4),ncol=count[i]*ndim,byrow=TRUE)
  sim_z_dist = array(NA,dim=c(niter,count[i],count[i]))
  sim_z[i,,1:(count[i]*ndim)] = temp
  for(iter in 1:niter){
    temp_mat = matrix(temp[iter,],count[i],ndim,byrow=TRUE)
    sim_z_dist[iter,,] = as.matrix(dist(temp_mat))
  }
  if(option == 1){
    pdf(oname3)
    for(rr in 2:count[i]){
      for(cc in 1:(rr-1)){
        plot(a,sim_z_dist[,rr,cc],type="l",main=paste("row:", rr, "col:", cc))
      }
    }
    dev.off()
  }
  
  # Draw the trace plots for pairwise distance of item latent positions for each school
  temp = matrix(scan(fname5),ncol=nitem*ndim,byrow=TRUE)
  sim_i_dist = array(NA,dim=c(niter,nitem,nitem))
  sim_i[i,,1:(nitem*ndim)] = temp
  for(iter in 1:niter){
    temp_mat = matrix(temp[iter,],nitem,ndim,byrow=TRUE)
    sim_i_dist[iter,,] = as.matrix(dist(temp_mat))
  }
  if(option == 1){
    pdf(oname4)
    for(rr in 2:nitem){
      for(cc in 1:(rr-1)){
        plot(a,sim_i_dist[,rr,cc],type="l",main=paste("row:", rr, "col:", cc))
      }
    }
    dev.off()
  }
  
  # Get $\sigma_{\epsilon_k}^2$ for each school
  temp = scan(fname3)
  sim_sigz[,i] = temp
}

#---------------------------------------------------------------------------------------------------------------
# Draw the trace plot $\sigma_{\epsilon_k}^2$ for each school
#---------------------------------------------------------------------------------------------------------------
pdf("TRACE/upper_sigma_z.pdf")
for(i in 1:nschool) plot(a,sim_sigz[,i],type="l",main=i)
dev.off();

#---------------------------------------------------------------------------------------------------------------
# Trace Plot for Hyperparameters
# "upper_gamma.pdf": Figure 5 of the supplemntary materials ($\gamma$). Hyperparameter for beta
# "upper_sigma_gamma.pdf": Figure 6 of the supplementary materials ($\sigma_{\beta}^2$). Variance of the hyperparameter for beta.
# "upper_sigma_d.pdf": Figure 7 of the supplementary material. Variances of random effects in the item dependence matrices #####
# "upper_mu.pdf": Figure 8 of the supplementary material. Fixed effects of the item dependence matrices.
# "upper_sigma_delta.pdf": Figure 9 of the supplementary material. Variance of the fixed effects of the item dependence matrices.
# "upper_school_dist.pdf": School distances calculated from the item dependence matrices
#---------------------------------------------------------------------------------------------------------------
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

sim_gamma  = matrix(scan("RESULT/sim_g1.log"),ncol=nitem*ncat,byrow=TRUE)
sim_varphi = matrix(scan("RESULT/sim_p1.log"),ncol=nitem*ncat,byrow=TRUE)
sim_delta  = matrix(scan("RESULT/sim_l1.log"),ncol=ndyad*ncat,byrow=TRUE)
sim_tau    = matrix(scan("RESULT/sim_u1.log"),ncol=ndyad*ncat,byrow=TRUE)
if(option==1){
  pdf("TRACE/upper_gamma.pdf")
  for(i in 1:nitem) plot(a, sim_gamma[,i], type="l", main=paste("tradition item ",i,sep=""))
  for(i in 1:nitem) plot(a, sim_gamma[,(nitem+i)], type="l", main=paste("innovation item ",i,sep=""))
  dev.off()
  pdf("TRACE/upper_sigma_beta.pdf")
  for(i in 1:nitem) plot(a, sim_varphi[,i], type="l", main=paste("tradition item ",i,sep=""))
  for(i in 1:nitem) plot(a, sim_varphi[,(nitem+i)], type="l", main=paste("innovation item ",i,sep=""))
  dev.off()
  pdf("TRACE/upper_mu.pdf")
  for(i in 1:ndyad) plot(a, sim_delta[,i], type="l", 
                         main=paste("tradition distance between item ",item_index[i,1]," and item ",item_index[i,2]))
  for(i in 1:ndyad) plot(a, sim_delta[,(ndyad+i)], type="l", 
                         main=paste("innovation distance between item ",item_index[i,1]," and item ",item_index[i,2]))
  dev.off()
  pdf("TRACE/upper_sigma_delta.pdf")
  for(i in 1:ndyad) plot(a, sim_tau[,i], type="l", 
                         main=paste("tradition distance between item ",item_index[i,1]," and item ",item_index[i,2]))
  for(i in 1:ndyad) plot(a, sim_tau[,(ndyad+i)], type="l", 
                         main=paste("innovation distance between item ",item_index[i,1]," and item ",item_index[i,2]))
  dev.off()
}

#---------------------------------------------------------------------------------------------------------------
# Reconstruction of the item latent positions using from $\mu$
# Figure 4 of the manuscript
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

pdf("PLOT/headmap_control.pdf")
heatmap(idist_control, symm=TRUE, keep.dendro=FALSE, cexRow=0.5, cexCol=0.5)
dev.off()

if(option == 1){
  ntrad = isoMDS(idist_control,k=2)
  ntrad_x = ntrad$points[,1]; ntrad_y = ntrad$points[,2];
  write.table(ntrad$points,"RSUM/position_control.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/position_control.txt"),ncol=2,byrow=TRUE)
  ntrad_x = temp[,1]; ntrad_y = temp[,2];
}

pdf("PLOT/mu_xy_tradition.pdf")
plot(ntrad_x, ntrad_y, xlab="", ylab="", type="n")
text(ntrad_x, ntrad_y, labels=1:nitem, cex=1.5)
dev.off()

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

pdf("PLOT/headmap_innovation.pdf");
heatmap(idist_innovation, symm=TRUE, keep.dendro=FALSE, cexRow=0.5, cexCol=0.5);
dev.off();

if(option == 1){
  ninno = isoMDS(idist_innovation,k=2)
  #ntrad = cmdscale(idist,eig=TRUE,k=2)
  ninno_x = ninno$points[,1]; ninno_y = ninno$points[,2];
  write.table(ninno$points,"RSUM/position_innovation.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/position_innovation.txt"),ncol=2,byrow=TRUE)
  ninno_x = temp[,1]; ninno_y = temp[,2];
}

pdf("PLOT/mu_xy_innovation.pdf")
plot(ninno_x, ninno_y, xlab="", ylab="", type="n")
text(ninno_x, ninno_y, labels=1:nitem, cex=1.5)
dev.off()

save.image("GEPS_REVERSE_MODEL2.RData")
