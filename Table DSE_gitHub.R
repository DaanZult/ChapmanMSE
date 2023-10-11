rm(list=ls())
gc()

set.seed(42)

library(matlib)
library(matrixStats)
library(brglm2)

# no. iterations
MC = 20000

results_scenarios = list()
scenarios = list()
scenarios[[1]] = list(N=100,   p_a=0.5,  p_b=0.2,  osn = 1)
scenarios[[2]] = list(N=100,   p_a=0.35, p_b=0.30, osn = 1)
scenarios[[3]] = list(N=500,   p_a=0.4,  p_b=0.15, osn = 1)
scenarios[[4]] = list(N=500,   p_a=0.25, p_b=0.20, osn = 1)
scenarios[[5]] = list(N=10000, p_a=0.3,  p_b=0.10, osn = 1)
scenarios[[6]] = list(N=10000, p_a=0.25, p_b=0.15, osn = 1)
scenarios[[7]] = list(N=100,   p_a=0.15, p_b=0.15, osn = 1)

generate_table = function(N = N, p_a = p_a, p_b=p_b, osn = osn){
  ### create dataset from contingency table
  list_a <- factor(c("in", "missed", "in", "missed"), levels=c("missed","in"))
  list_b <- factor(c("in", "missed", "missed", "in"), levels=c("missed","in"))
  count <- c((p_a*p_b*N), ((1-p_a)*(1-p_b)*N),
            ((p_a)*(1-p_b)*N),((1-p_a)*(p_b)*N))
  
  ### dataset of frequencies
  data_1 <- data.frame(list_a, list_b, count)
  ### create offset variable
  data_1$offset <- c(log(osn), log(1), log(1), log(1))
  ### model estimated counts from independence model and include offset term
  GLM <- glm(formula = count ~ list_a + list_b,
                  family = poisson(link = "log"),
                  data = data_1,
                  offset = data_1$offset)
  
  summary(GLM)
  GLM$fitted
  data_1$est_countoff<- c(GLM$fitted)
  ### check to see if offset gave or = osn
  or_offset<- (GLM$fitted[1]*GLM$fitted[2])/
    (GLM$fitted[3]*GLM$fitted[4])
  or_offset
  ### check the response of the lists again
  ### reponse should be p_a
  list_a_resp <- (GLM$fitted[1]+GLM$fitted[3])/
    (GLM$fitted[1]+GLM$fitted[2]+
       GLM$fitted[3]+GLM$fitted[4])
  list_a_resp
  
  ### response should be p_b
  list_b_resp <- (GLM$fitted[1]+GLM$fitted[4])/
    (GLM$fitted[1]+GLM$fitted[2]+
       GLM$fitted[3]+GLM$fitted[4])
  list_b_resp
  ### total size
  data_1$total_samp <- sum(data_1$est_countoff)
  data_1_final <- cbind(data_1,data_1$total_samp)
  
  
  data_1_final$mysum <- rmultinom(1, N, data_1_final$est_countoff)
  data_1_final$mysum <- as.numeric(data_1_final$mysum)
  table <- data_1_final[!(data_1_final$list_a=="missed"
                         &data_1_final$list_b=="missed"),]
  return(table)
}

estimate_DSEs = function(table = table){
  
  table_save = table
  n=sum(table$mysum)
  ### model estimated counts from independence model
  GLM <- glm(formula = mysum ~ list_a + list_b,
                       family = poisson(link = "log"),
                       data = table)
  
  m00_LP <- exp(summary(GLM)$coefficients[1,1])
  N_LP = n + m00_LP
  
  # Bailey adjustment
  table[table$list_a=="in"&table$list_b=="in","mysum"]     = table[table$list_a=="in"&table$list_b=="in","mysum"]+1
  table[table$list_a=="missed"&table$list_b=="in","mysum"] = table[table$list_a=="missed"&table$list_b=="in","mysum"]-1
  
  ### model estimated counts from independence model
  GLM <- glm(formula = mysum ~ list_a + list_b,
                  family = poisson(link = "log"),
                  data = table)
  
  m00_Bailey <- exp(summary(GLM)$coefficients[1,1])
  N_Bailey = n + m00_Bailey

  
  table = table_save
  
  # Evans adjustment
  table[,"mysum"] = table[,"mysum"]+0.5
  
  GLM <- glm(formula = mysum ~ list_a + list_b,
                       family = poisson(link = "log"),
                       data = table)
  
  m00_Evans <- exp(summary(GLM)$coefficients[1,1])
  N_Evans = n + m00_Evans
  
  table = table_save
  
  GLM <- glm(formula = mysum ~ list_a + list_b,
                  family = poisson(link = "log"),
                  data = table, method = "brglmFit",type="MPL_Jeffreys")
  
  m00_Firth <- exp(summary(GLM)$coefficients[1,1])
  N_Firth = n + m00_Firth
  
  GLM <- glm(formula = mysum ~ list_a + list_b,
                  family = poisson(link = "log"),
                  data = table, method = "brglmFit", type = "AS_mean")
  
  m00_Kosmidis <- exp(summary(GLM)$coefficients[1,1])
  N_Kosmidis = n + m00_Kosmidis
  
  GLM <- try(glm(formula = mysum ~ list_a + list_b,
                      family = poisson(link = "log"),
                      data = table, method = "brglmFit", type = "correction"), silent = TRUE)
  summary(GLM)
  if (class(GLM)[1]=="try-error"){
    N_Cordeiro = N_Firth}
  if (class(GLM)[1]!="try-error"){
    m00_Cordeiro <- exp(summary(GLM)$coefficients[1,1])
    N_Cordeiro = n + m00_Cordeiro
  }
  
  table = table_save
  
  # Chapman adjustment
  table[table$list_a=="in"&table$list_b=="in","mysum"] = table[table$list_a=="in"&table$list_b=="in","mysum"]+1
  
  GLM <- glm(formula = mysum ~ list_a + list_b,
                       family = poisson(link = "log"),
                       data = table)
  
  m00_Chapman <- exp(summary(GLM)$coefficients[1,1])
  N_Chapman = n + m00_Chapman
  
  
  return(list(n, N_LP, N_Bailey, N_Evans,
              N_Firth, N_Cordeiro, N_Kosmidis,
              N_Chapman))
}

### to store estimates
for (s in 1:length(scenarios)){
p_a = scenarios[[s]]$p_a
p_b = scenarios[[s]]$p_b
N = scenarios[[s]]$N
osn <- 1
n = NULL
N_LP <- NULL
N_Bailey <- NULL
N_Evans <- NULL

N_Firth <- NULL
N_Cordeiro <- NULL
N_Kosmidis <- NULL

N_Chap <- NULL

for (mc in 1:MC){
  table = generate_table(N = N, p_a = p_a, p_b = p_b, osn = osn)
  DSEs = estimate_DSEs(table)
  
  ### store estimates in total_est
  n = c(n, DSEs[[1]])
  
  N_LP       <- c(N_LP,       DSEs[[2]])
  N_Bailey   <- c(N_Bailey,   DSEs[[3]])
  N_Evans    <- c(N_Evans,    DSEs[[4]])
  
  N_Firth    <- c(N_Firth,    DSEs[[5]])
  N_Cordeiro <- c(N_Cordeiro, DSEs[[6]])
  N_Kosmidis <- c(N_Kosmidis, DSEs[[7]])
  
  N_Chap     <- c(N_Chap,DSEs[[8]])
}

results_scenarios[[s]]=as.matrix(cbind(n,
                                       N_LP,N_Bailey,N_Evans,
                                       N_Firth, N_Cordeiro, N_Kosmidis,
                                       N_Chap))

print(s)
print(round(colMeans(results_scenarios[[s]]),1))
}

for (s in 1:length(scenarios)){
  print(sum(results_scenarios[[s]][,"N_LP"]>(10*mean(results_scenarios[[s]][,"N_Chap"]))))
  results_scenarios[[s]][results_scenarios[[s]][,"N_LP"]>10*mean(results_scenarios[[s]][,"N_Chap"]),"N_LP"] = 
    max(results_scenarios[[s]][,"N_Chap"])
}

options(scipen=99)

est = numeric()
est_sds = numeric()
est_pvals = numeric()
for (s in 1:length(scenarios)){
  est_sd = cbind(t(t(format(round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1),nsmall=1))),t(t(round(colSds(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE)/MC^0.5,2))))
  est_sd = apply(est_sd,1,paste0, collapse = " (")
  est_sd = paste0(est_sd,")")
  est_sds = cbind(est_sds,est_sd)
  TTESTS = apply(results_scenarios[[s]][,c(1,2,3,4,8)]-scenarios[[s]]$N,2,t.test)
  est_pval = cbind(t(t(format(round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1),nsmall=1))),as.numeric(unlist(TTESTS)[c(3,14,25,36,47)]))
  est_pval[as.numeric(est_pval[,2])<=0.001,2] = "***"
  est_pval[as.numeric(est_pval[,2])>0.001&est_pval[,2]<=0.01,2] = "**"
  est_pval[as.numeric(est_pval[,2])>0.01&est_pval[,2]<=0.05,2] = "*"
  est_pval[as.numeric(est_pval[,2])>0.05,2] = ""
  est_pval = apply(est_pval,1,paste0, collapse = "")
  est_pvals = cbind(est_pvals,est_pval)
  est = cbind(est,t(t(round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1))))
}

X = cbind(1:7,
          c(100,100,500,500,10000,10000,100),
          c(0.5,0.35,0.4,0.25,0.3,0.25,0.15),
          c(0.2,0.3,0.15,0.2,0.1,0.15,0.15),
          t(est_pvals))
DSE_means = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(DSE_means)
DSE_means = gsub(".","&.",DSE_means,fixed=TRUE)
DSE_means = gsub("  ","",DSE_means)
print(DSE_means)

est_sds = numeric()
for (s in 1:length(scenarios)){
  est_sd  = format(round(colSds(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1),nsmall=1)
  est_sds = cbind(est_sds, est_sd)
}

est_sds = est_sds[-1,]

X = cbind(1:7,
          c(100,100,500,500,10000,10000,100),
          c(0.5,0.35,0.4,0.25,0.3,0.25,0.15),
          c(0.2,0.3,0.15,0.2,0.1,0.15,0.15),
          t(est_sds))
DSE_sds = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(DSE_sds)
DSE_sds = gsub(".","&.",DSE_sds,fixed=TRUE)
DSE_sds = gsub("  ","",DSE_sds)
print(DSE_sds)

est_rmses = numeric()
for (s in 1:length(scenarios)){
  est_rmse  = format(round(colMeans(((results_scenarios[[s]][,c(1,2,3,4,8)]-scenarios[[s]]$N)^2),na.rm=TRUE)^0.5,1),nsmall=1)
  est_rmses = cbind(est_rmses, est_rmse)
}

est_rmses = est_rmses[-1,]

X = cbind(1:7,
          c(100,100,500,500,10000,10000,100),
          c(0.5,0.35,0.4,0.25,0.3,0.25,0.15),
          c(0.2,0.3,0.15,0.2,0.1,0.15,0.15),
          t(est_rmses))
DSE_rmses = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(DSE_rmses)
DSE_rmses = gsub(".","&.",DSE_rmses,fixed=TRUE)
DSE_rmses = gsub("  ","",DSE_rmses)
print(DSE_rmses)

est_sds = numeric()
est_rmses = numeric()
for (s in 1:length(scenarios)){
  est_sd  = format(round(colSds(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1),nsmall=1)
  est_sds = cbind(est_sds, est_sd)
}

for (s in 1:length(scenarios)){
  est_rmse  = format(round(colMeans(((results_scenarios[[s]][,c(1,2,3,4,8)]-scenarios[[s]]$N)^2),na.rm=TRUE)^0.5,1),nsmall=1)
  est_rmses = cbind(est_rmses, est_rmse)
}

est_sds_rmses = rbind(est_sds[-1,],est_rmses[-1,])

X = cbind(1:7,
          t(est_sds_rmses))
DSE_sds_rmses = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(DSE_sds_rmses)
DSE_sds_rmses = gsub(".","&.",DSE_sds_rmses,fixed=TRUE)
DSE_sds_rmses = gsub("  ","",DSE_sds_rmses)
print(DSE_sds_rmses)

DSE_means
DSE_sds
DSE_rmses
DSE_sds_rmses
