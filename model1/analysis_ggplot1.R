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
library(gplots)
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

sampvec = rep(NA,sum(count)); index = 0;
for(i in 1:nschool){
  sampvec[(index+1):(index+count[i])] = sampsum[i,(1:count[i])]
  index = index + count[i]
}
table(sampvec)

#---------------------------------------------------------------------------------------------------------------
# Proportion of positive response of innnovation and regular.
# Figure 1 of the manuscript
#---------------------------------------------------------------------------------------------------------------

renov <- read.table("DATA/renov.txt")
temp_sum <- data.frame(itemsum, renov, count)
# proportion
temp_sum2 <- data.frame(temp_sum %>% group_by(V1) %>% summarise_all(sum))
# reshape
temp_sum3 <- data.frame(temp_sum2[, 1], temp_sum2[,2:73] / temp_sum2[,74])
colnames(temp_sum3) <- c("renov", 1:72)
temp_sum4 <- melt(data = temp_sum3, id.vars = "renov")
pdf("GGPLOT/prop_positive.pdf")
gg <- ggplot(temp_sum4, aes(x = variable, y = value * 100)) +
  geom_point(aes(color = as.factor(renov), shape = as.factor(renov)), size = 2.5) +
  scale_x_discrete(breaks = round(seq(from = 0, to = 72, by = 10))) + 
  xlab("Items") + ylab("Proportion or Positive Responses") +
  scale_shape_discrete(labels=c("Regular", "Innovation")) +
  scale_colour_discrete(labels=c("Regular", "Innovation")) +
  theme(legend.title = element_blank(), legend.position = "top")
print(gg)
dev.off()

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

#---------------------------------------------------------------------------------------------------------------
# Figure 2 of the Supplementary Materials 
# (Three number summary of the person intercept parameter estimate categorized by the total scores)
#---------------------------------------------------------------------------------------------------------------
pdf("PLOT/hpd_theta.pdf")
sampvec = rep(NA,sum(count)); index = 0;
for(i in 1:nschool){
  sampvec[(index+1):(index+count[i])] = sampsum[i,(1:count[i])]
  index = index + count[i]
}

theta_int = 10:66
theta_mean = matrix(NA,nitem+1,5)
for(i in theta_int){
  temp = theta_result[sampvec==i]
  theta_mean[(i+1),] = quantile(temp,c(0.025,0.25,0.5,0.75,0.975))
}
x.axis = theta_int
y.axis = x.axis

plot(x.axis,y.axis,ylim=c(-4.0,4.0),xlab="total score",ylab=expression(theta),type="n")
for(i in theta_int){
  points(i,theta_mean[i,3],pch=20,col=4)
  lines(c(i,i),c(theta_mean[i,1],theta_mean[i,5]),lwd=2)
} 
for(i in 1:5) abline(v=i*10, lty=3, col=2)
dev.off()

pdf("GGPLOT/hpd_theta.pdf")
temp <- data.frame(theta_mean, total_score = 1:73)[11:66,]
gg <- ggplot(temp)+
  geom_linerange(aes(x=total_score, ymin=X1, ymax=X5), size = 5) +
  geom_point(aes(x=total_score, y=X3), color = "blue", size = 0.5) +
  geom_vline(xintercept=round(seq(from = 10, to = 66, by = 10)), linetype = 'dotted', color='red', size = 0.7) +
  scale_x_continuous(breaks = round(seq(from = 10, to = 66, by = 10))) +
  xlab("Total score") + ylab(expression(theta))
print(gg)
dev.off()

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
# Histogram of proportion of positive response of each school.
# Figure 8 (a, b) of the manuscript 
# Draw for all schools
#---------------------------------------------------------------------------------------------------------------

temp_hist <- sampsum/72
pdf(paste("GGPLOT/prop_agree.pdf", sep = ""))
for(i in 1:72){
  temp_prob2 <- data.frame(value = temp_hist[i, ])
  temp_prob2 <- data.frame(value = temp_prob2[is.na(temp_prob2)!=1, ])
  hist_temp <- ggplot(temp_prob2, aes(x=value)) + 
    theme_minimal()+ geom_histogram(binwidth=0.05,color="black", fill="white") +
    xlim(c(0, 1)) +
    xlab("Proportion of Agree") +
    ylab("Frequency")
  print(hist_temp)
}
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Proportion of positive response of each school.
# Figure 8 (c) of the manuscript
#---------------------------------------------------------------------------------------------------------------

temp_prob <- as.data.frame(itemprob[c(10, 28), ])
colnames(temp_prob) <- c(1:72)
temp_prob$school <- c(10, 28)
temp_prob2 <- melt(data = temp_prob, id.vars = "school")
pdf("GGPLOT/prop_positive_each.pdf")
gg <- ggplot(temp_prob2, aes(x = variable, y = value * 100, group = school)) +
  geom_line(aes(linetype = as.factor(school), color = as.factor(school)))+
  geom_point(aes(color = as.factor(school), shape = as.factor(school))) +
  scale_x_discrete(breaks = round(seq(from = 0, to = 72, by = 10))) + 
  xlab("Items") + ylab("Proportion or Positive Responses") +
  # scale_shape_discrete(labels=c("10", "28")) +
  scale_colour_discrete(labels=c("10", "28")) +
  theme(legend.title = element_blank(), legend.position = "top")
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Mean and 95% HPD intervals for $\gamma$
# Figure 2(a) of the manuscript
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

pdf("GGPLOT/hpd_gamma.pdf")
gamma_result <- data.frame(gamma_result)
gamma_result$item <- 1:72
ggplot(gamma_result)+
  geom_linerange(aes(x=item, ymin=X2, ymax=X3), size = 1) +
  geom_point(aes(x=item, y=X1), color = "blue") +
  geom_vline(xintercept=round(seq(from = 0, to = 72, by = 10)), linetype = 'dotted', color='red', size = 0.7) +
  scale_x_continuous(breaks = round(seq(from = 0, to = 72, by = 10))) +
  xlab("Item") + ylab(expression(gamma))
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Mean and 95% HPD intervals for $\sigma_{\beta}^2$
# Figure 1 of the supplementary materials
#---------------------------------------------------------------------------------------------------------------
varphi_mean = apply(sim_varphi,2,mean)
varphi_mcmc = mcmc(sim_varphi)
varphi_result = matrix(NA,length(varphi_mean),3)
varphi_result[,1] = varphi_mean
varphi_result[,2:3] = HPDinterval(varphi_mcmc,prob=0.95)
round(varphi_result,4)
min(varphi_result)
max(varphi_result)

pdf("PLOT/hpd_sigma_beta.pdf")
x.axis = c(1:nitem)
y.axis = x.axis
plot(x.axis,y.axis,ylim=c(0.25,2.00),xlab="item",ylab=expression(sigma[beta]^2),type="n")
for(i in 1:nitem){
  lines(c(i,i),varphi_result[i,2:3],lwd=2)
  points(i,varphi_result[i,1],pch=20,col=4)
} 
for(i in 1:7) abline(v=i*10, lty=3, col=2)
dev.off()

pdf("GGPLOT/hpd_sigma_beta.pdf")
varphi_result <- data.frame(varphi_result)
varphi_result$item <- 1:72
gg <- ggplot(varphi_result)+
  geom_linerange(aes(x=item, ymin=X2, ymax=X3), size = 1) +
  geom_point(aes(x=item, y=X1), color = "blue") +
  geom_vline(xintercept=round(seq(from = 0, to = 72, by = 10)), linetype = 'dotted', color='red', size = 0.7) +
  scale_x_continuous(breaks = round(seq(from = 0, to = 72, by = 10))) +
  xlab("Item") + ylab(expression(sigma[beta]^2))
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Mean and 95% HPD intervals for $\sigma_z^2$.
# Figure 3 of the supplementary materials
#---------------------------------------------------------------------------------------------------------------
sigz_mean = apply(sim_sigz,2,mean)
sigz_mcmc = mcmc(sim_sigz)
sigz_result = matrix(NA,length(sigz_mean),3)
sigz_result[,1] = sigz_mean
sigz_result[,2:3] = HPDinterval(sigz_mcmc,prob=0.95)
round(sigz_result,4)
min(sigz_result)
max(sigz_result)

pdf("PLOT/hpd_sigma_z.pdf")
x.axis = c(1:nschool)
y.axis = x.axis
plot(x.axis,y.axis,ylim=c(0.0,15.0),xlab="school",ylab=expression(sigma[z]^2),type="n")
for(i in 1:nschool){
  lines(c(i,i),sigz_result[i,2:3],lwd=2)
  points(i,sigz_result[i,1],pch=20,col=4)
} 
for(i in 1:6) abline(v=i*10, lty=3, col=2)
dev.off()

# Correction
pdf("GGPLOT/hpd_sigma_z.pdf")
sigz_result <- data.frame(sigz_result)
sigz_result$school <- 1:nschool
gg <- ggplot(sigz_result)+
  geom_linerange(aes(x=school, ymin=X2, ymax=X3), size = 1) +
  geom_point(aes(x=school, y=X1), color = "blue") +
  geom_vline(xintercept=round(seq(from = 0, to = nschool, by = 10)), linetype = 'dotted', color='red', size = 0.7) +
  scale_x_continuous(breaks = round(seq(from = 0, to = nschool, by = 10))) +
  xlab("School") + ylab(expression(sigma[z]^2))
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Visualization of the item latent positions using from estimated item distance matrix at the upper level.
# 1. Heatmap for estimated $\mu$ (Figure 4 of the Manuscript)
# 2. Reconstruction of the item latent positions using from $\mu$ (Figure 3 of the Manuscript)
#---------------------------------------------------------------------------------------------------------------
temp = matrix(scan("RESULT/sum_l1.log"),ncol=ndyad,byrow=TRUE)
temp1 = temp[1,1:ndyad];
idist = matrix(0,nitem,nitem);
logidist = matrix(0,nitem,nitem);
index = 0;
for(i in 2:nitem){
  for(j in 1:(i-1)){
    index = index + 1
    idist[i,j] = exp(temp1[index])
    logidist[i,j] = temp1[index]
  }
}
idist = idist + t(idist)
logidist = logidist + t(logidist)

#---------------------------------------------------------------------------------------------------------------
# Heatmap for estimated $\mu$ 
# Figure 4 of the Manuscript
#---------------------------------------------------------------------------------------------------------------
pdf("GGPLOT/heatmap_mu.pdf");
heatmap(exp(-idist), symm=TRUE, keep.dendro=FALSE)
dev.off();

#---------------------------------------------------------------------------------------------------------------
# Reconstruction of the item latent positions using from $\mu$
# Figure 3 of the manuscript
#---------------------------------------------------------------------------------------------------------------
if(option == 1){
  ntrad = isoMDS(idist,k=2)
  ntrad_x = ntrad$points[,1]; ntrad_y = ntrad$points[,2];
  write.table(ntrad$points,"RSUM/position.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/position.txt"),ncol=2,byrow=TRUE)
  ntrad_x = temp[,1]; ntrad_y = temp[,2];
}

pdf("PLOT/position_mu.pdf")
plot(ntrad_x, ntrad_y, xlab="", ylab="", type="n")
text(ntrad_x, ntrad_y, labels=1:nitem, cex=1.5)
dev.off()

pdf("GGPLOT/position_mu.pdf")
temp_ls <- data.frame(x = ntrad_x, y = ntrad_y)
gg <- ggplot() + 
  geom_text(data = temp_ls, aes(x = ntrad_x, y = ntrad_y, label=1:nrow(temp_ls)), size = 5, fontface = "bold") +
  xlab(NULL) + ylab(NULL)
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Visualization of the school latent positions using from estimated item distance matrix at the upper level.
# Figure 5 of the manuscript
#---------------------------------------------------------------------------------------------------------------
sdist = matrix(scan("RESULT/sum_a1.log"),ncol=nschool,byrow=TRUE)
if(option == 1){
  spoint = isoMDS(sdist, k=2)
  spoint_x = spoint$points[,1]; spoint_y = spoint$points[,2];
  write.table(spoint$points,"RSUM/school_pos.txt",row.names=FALSE,col.names=FALSE)
}else{
  temp = matrix(scan("RSUM/school_pos.txt"),ncol=2,byrow=TRUE)
  spoint_x = temp[,1]; spoint_y = temp[,2];
}

pdf("PLOT/position_school.pdf")
plot(spoint_x, spoint_y, type="n")
text(spoint_x, spoint_y, labels=1:nschool, cex=1.5)
dev.off()

pdf("PLOT/position_school_renov.pdf")
renov = scan("DATA/renov.txt")
plot(spoint_x, spoint_y, xlab="", ylab="", type="n")
text(spoint_x, spoint_y, labels=1:nschool, col=renov+1, cex=1.5)
dev.off()
length(renov)

pdf("GGPLOT/position_school.pdf")
temp_ls <- data.frame(x = spoint_x, y = spoint_y )
gg <- ggplot() + 
  geom_text(data = temp_ls, aes(x = spoint_x, y = spoint_y, label=1:nrow(temp_ls)) , size = 5, fontface = "bold") +
  xlab(NULL) + ylab(NULL)
print(gg)
dev.off()

pdf("GGPLOT/position_school_renov.pdf")
renov = read.table("DATA/renov.txt")
temp_ls <- data.frame(x = spoint_x, y = spoint_y , z = renov)
gg <- ggplot() + 
  geom_text(data = temp_ls, aes(x = spoint_x, y = spoint_y, label=1:nrow(temp_ls), colour = factor(V1)) , size = 5, fontface = "bold")  +
  theme(legend.title = element_blank(), legend.position = "top")+
  xlab(NULL) + ylab(NULL)+
  scale_shape_discrete(labels=c("Regular", "Innovation")) +
  scale_colour_discrete(labels=c("Regular", "Innovation")) +
  theme(legend.title = element_blank())
print(gg)
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Item Latent Spaces for Each School
# Figure 6 of the manuscript / Figure 9-24 of the supplementary material
#---------------------------------------------------------------------------------------------------------------
pdf("PLOT/item_ls.pdf")
item_ls = array(0, dim=c(nschool, nitem, ndim))
item_recover = array(0, dim=c(nschool, nitem, ndim))
for(k in 1:nschool){
  if(k < 10) fname1 = paste("RESULT/sum0",k,"_i1.log",sep="")
  else fname1 = paste("RESULT/sum",k,"_i1.log",sep="")
  temp = as.matrix(read.table(fname1))
  item_ls[k,,] = matrix(temp[1,],ncol=ndim,byrow=TRUE)
  plot(item_ls[k,,], xlab="", ylab="", type="n")
  text(item_ls[k,,], labels = 1:nitem)
  print(k)
}
dev.off()


samp_ls <- list()
for(k in 1:nschool){
  if(k < 10) fname1 = paste("RESULT/sum0",k,"_z1.log",sep="")
  else fname1 = paste("RESULT/sum",k,"_z1.log",sep="")
  temp = as.matrix(read.table(fname1))
  samp_ls[[k]] = matrix(temp[1,],ncol=ndim,byrow=TRUE)
  # plot(samp_ls[[k]][,1] ~ samp_ls[[k]][,2], xlab="", ylab="")
  print(k)
}

library("ggrepel")
pdf("GGPLOT/item_ls.pdf")
for(k in 1:nschool){
  temp_ls <- data.frame(item_ls[k,,])
  temp_ls2 <- data.frame(samp_ls[[k]])
  gg <- ggplot() +
    geom_point(data = temp_ls2,  aes(x = X1, y = X2), col = 1) +
    geom_text(data = temp_ls, aes(x = X1, y = X2, label=1:nrow(temp_ls)), size = 4.5, fontface = "bold", col = 2) +
    xlab(NULL) + ylab(NULL)
  print(gg)
}
dev.off()

#---------------------------------------------------------------------------------------------------------------
# Model Assessment of Y and U for Each School
# Table 10 and 11 of the supplementary material
#---------------------------------------------------------------------------------------------------------------
auc.Y.final = rep(NA,nschool); sen.Y.final = rep(NA,nschool) 
spe.Y.final = rep(NA,nschool); ove.Y.final = rep(NA,nschool)
auc.U.final = rep(NA,nschool); sen.U.final = rep(NA,nschool)
spe.U.final = rep(NA,nschool); ove.U.final = rep(NA,nschool)
phi.Y.final = rep(NA,nschool); phi.U.final = rep(NA,nschool)

TN.Y.sum = rep(0,nschool); TN.U.sum = rep(0,nschool);
TP.Y.sum = rep(0,nschool); TP.U.sum = rep(0,nschool);
FN.Y.sum = rep(0,nschool); FN.U.sum = rep(0,nschool);
FP.Y.sum = rep(0,nschool); FP.U.sum = rep(0,nschool);

w_temp = array(0, dim=c(nschool, nitem, nitem))
for(k in 1:nschool){
  if(k < 10) fname1 = paste("RESULT/sum0",k,"_d1.log",sep="")
  else fname1 = paste("RESULT/sum",k,"_d1.log",sep="")
  temp = as.matrix(read.table(fname1))
  index = 0
  for(i in 2:nitem){
    for(j in 1:(i-1)){
      index = index + 1
      w_temp[k,i,j] = temp[1,index]
    }
  }
  w_temp[k,,] = w_temp[k,,] + t(w_temp[k,,])
}

for(a in 1:nschool){
  beta_temp = apply(sim_beta[a,,], 2, mean)
  theta_temp = apply(sim_theta[a,,1:count[a]], 2, mean)
  z_temp = array(0, dim=c(niter,count[a],ndim))
  z_dist_temp = array(0, dim=c(niter,count[a],count[a]))
  z_dist = matrix(0,count[a],count[a])
  for(iter in 1:niter){
    z_temp[iter,,] = matrix(sim_z[a,iter,1:(count[a]*ndim)], ncol=ndim, byrow=TRUE)
    z_dist_temp[iter,,] = as.matrix(dist(z_temp[iter,,]))
  }
  for(i in 1:count[a]){
    for(j in 1:count[a]){
      z_dist[i,j] = mean(z_dist_temp[,i,j])
    }
  }
  
  lower_Y_orig = {}
  lower_Y_pred = {}
  for(k in 1:nitem){
    temp_mat_Y = Y[a,k,1:count[a],1:count[a]]
    lower_Y_orig = c(lower_Y_orig, temp_mat_Y[lower.tri(temp_mat_Y)])
    prob_mat_Y_pred = plogis(beta_temp[k] - z_dist)
    temp_mat_Y_pred = matrix(rbinom(count[a]*count[a],1,prob_mat_Y_pred),count[a],count[a])
    lower_Y_pred = c(lower_Y_pred, temp_mat_Y_pred[lower.tri(temp_mat_Y_pred)])
  }
  
  lower_U_orig = {}
  lower_U_pred = {}
  for(k in 1:count[a]){
    temp_mat_U = U[a,k,,]
    lower_U_orig = c(lower_U_orig, temp_mat_U[lower.tri(temp_mat_U)])
    prob_mat_U_pred = plogis(theta_temp[k] - w_temp[a,,])
    temp_mat_U_pred = matrix(rbinom(nitem*nitem,1,prob_mat_U_pred),nitem,nitem)
    lower_U_pred = c(lower_U_pred, temp_mat_U_pred[lower.tri(temp_mat_U_pred)])
  }
  
  temp.table = table(predicted = lower_Y_pred, actual = lower_Y_orig)
  TN.Y.final = temp.table[1,1] / sum(temp.table)
  TP.Y.final = temp.table[2,2] / sum(temp.table)
  FP.Y.final = temp.table[2,1] / sum(temp.table)
  FN.Y.final = temp.table[1,2] / sum(temp.table)
  TN.Y.sum[a] = TN.Y.sum[a] + TN.Y.final
  TP.Y.sum[a] = TP.Y.sum[a] + TP.Y.final
  FN.Y.sum[a] = FN.Y.sum[a] + FN.Y.final
  FP.Y.sum[a] = FP.Y.sum[a] + FP.Y.final
  spe.Y.final[a] = TN.Y.final /(TN.Y.final+FP.Y.final)
  sen.Y.final[a] = TP.Y.final /(FN.Y.final+TP.Y.final)
  ove.Y.final[a] = (TP.Y.final + TN.Y.final)
  phi.Y.final[a] = (TP.Y.final * TN.Y.final + FP.Y.final * FN.Y.final) / sqrt((TP.Y.final + FN.Y.final) * (TN.Y.final + FP.Y.final) * (TP.Y.final + FP.Y.final) * (TN.Y.final + FN.Y.final))
  temp = roc(lower_Y_orig, lower_Y_pred)
  auc.Y.final[a] = auc(temp)
  
  temp.table = table(predicted = lower_U_pred, actual = lower_U_orig)
  TN.U.final = temp.table[1,1] / sum(temp.table)
  TP.U.final = temp.table[2,2] / sum(temp.table)
  FP.U.final = temp.table[2,1] / sum(temp.table)
  FN.U.final = temp.table[1,2] / sum(temp.table)
  TN.U.sum[a] = TN.U.sum[a] + TN.U.final
  TP.U.sum[a] = TP.U.sum[a] + TP.U.final
  FN.U.sum[a] = FN.U.sum[a] + FN.U.final
  FP.U.sum[a] = FP.U.sum[a] + FP.U.final
  spe.U.final[a] = TN.U.final /(TN.U.final+FP.U.final)
  sen.U.final[a] = TP.U.final /(FN.U.final+TP.U.final)
  ove.U.final[a] = (TP.U.final + TN.U.final)
  phi.U.final[a] = (TP.U.final * TN.U.final + FP.U.final * FN.U.final) / sqrt((TP.U.final + FN.U.final) * (TN.U.final + FP.U.final) * (TP.U.final + FP.U.final) * (TN.U.final + FN.U.final))
  temp = roc(lower_U_orig, lower_U_pred)
  auc.U.final[a] = auc(temp)
  
  print(c(a,round(spe.Y.final[a],3),round(spe.U.final[a],3),round(sen.Y.final[a],3),round(sen.U.final[a],3),round(phi.Y.final[a],3),round(phi.U.final[a],3),round(auc.Y.final[a],3),round(auc.U.final[a],3)))
}

final.Y = cbind(sen.Y.final,spe.Y.final,ove.Y.final,phi.Y.final,auc.Y.final)
final.U = cbind(sen.U.final,spe.U.final,ove.U.final,phi.U.final,auc.U.final)
final.Y.mean = apply(final.Y,2,mean)
final.U.mean = apply(final.U,2,mean)
final.mean = rbind(final.Y.mean,final.U.mean)
colnames(final.Y) = c("Sensitivity", "Specificity", "Overall", "MCC", "AUC")
colnames(final.U) = c("Sensitivity", "Specificity", "Overall", "MCC", "AUC")
colnames(final.mean) = c("Sensitivity", "Specificity", "Overall", "MCC", "AUC")
write.table(round(final.Y,3),"RSUM/final_Y.txt",row.names=FALSE)
write.table(round(final.U,3),"RSUM/final_U.txt",row.names=FALSE)
write.table(round(final.mean,3),"RSUM/final_mean.txt",row.names=FALSE)

#---------------------------------------------------------------------------------------------------------------
# Save Computation Result
#---------------------------------------------------------------------------------------------------------------
save.image("GEPS_REVERSE_MODEL1.RData")
