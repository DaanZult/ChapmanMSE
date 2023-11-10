rm(list=ls())
gc()

set.seed(42)

library(matlib)
library(matrixStats)
library(brglm2)

# no. iterations
MC = 20000

# Define scenario's
results_scenarios = list()
scenarios = list()
scenarios[[1]] = list(N=100,   p_a=0.5,  p_b=0.2,  osn = 1)
scenarios[[2]] = list(N=100,   p_a=0.35, p_b=0.30, osn = 1)
scenarios[[3]] = list(N=500,   p_a=0.4,  p_b=0.15, osn = 1)
scenarios[[4]] = list(N=500,   p_a=0.25, p_b=0.20, osn = 1)
scenarios[[5]] = list(N=10000, p_a=0.3,  p_b=0.10, osn = 1)
scenarios[[6]] = list(N=10000, p_a=0.25, p_b=0.15, osn = 1)
scenarios[[7]] = list(N=100,   p_a=0.15, p_b=0.15, osn = 1)

# Define function to generate contingency table 
# (see Hammond et al., 2023 for further details)
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

# Function to obtain different population size estimates

estimate_DSEs = function(table = table){
  
  table_save = table
  n=sum(table$mysum)
  ### model estimated counts from independence model
  GLM <- glm(formula = mysum ~ list_a + list_b,
             family = poisson(link = "log"),
             data = table)
  
  m00_LP <- exp(summary(GLM)$coefficients[1,1])
  N_LP = n + m00_LP
  
  # Obtain Bailey estimate
  table[table$list_a=="in"&table$list_b=="in","mysum"]     = table[table$list_a=="in"&table$list_b=="in","mysum"]+1
  table[table$list_a=="missed"&table$list_b=="in","mysum"] = table[table$list_a=="missed"&table$list_b=="in","mysum"]-1
  
  ### model estimated counts from independence model
  GLM <- glm(formula = mysum ~ list_a + list_b,
             family = poisson(link = "log"),
             data = table)
  
  m00_Bailey <- exp(summary(GLM)$coefficients[1,1])
  N_Bailey = n + m00_Bailey
  
  
  table = table_save
  
  # Obtain Evans estimate 
  table[,"mysum"] = table[,"mysum"]+0.5
  
  GLM <- glm(formula = mysum ~ list_a + list_b,
             family = poisson(link = "log"),
             data = table)
  
  m00_Evans <- exp(summary(GLM)$coefficients[1,1])
  N_Evans = n + m00_Evans
  
  table = table_save
  
  # Obtain Firth estimate
  GLM <- glm(formula = mysum ~ list_a + list_b,
             family = poisson(link = "log"),
             data = table, method = "brglmFit",type="MPL_Jeffreys")
  
  m00_Firth <- exp(summary(GLM)$coefficients[1,1])
  N_Firth = n + m00_Firth
  
  # Obtain Kosmidis estimate
  GLM <- glm(formula = mysum ~ list_a + list_b,
             family = poisson(link = "log"),
             data = table, method = "brglmFit", type = "AS_mean")
  
  m00_Kosmidis <- exp(summary(GLM)$coefficients[1,1])
  N_Kosmidis = n + m00_Kosmidis
  
  # Obtain Cordeiro estimate
  GLM <- try(glm(formula = mysum ~ list_a + list_b,
                 family = poisson(link = "log"),
                 data = table, method = "brglmFit", type = "correction"), silent = TRUE)
  if (class(GLM)[1]=="try-error"){
    N_Cordeiro = N_Firth}
  if (class(GLM)[1]!="try-error"){
    m00_Cordeiro <- exp(summary(GLM)$coefficients[1,1])
    N_Cordeiro = n + m00_Cordeiro
  }
  
  table = table_save
  
  # Obtain Chapman estimate
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

### Loop over scenario 1 to 7
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

# Start Monte Carlo simulation for scenario s  
  for (mc in 1:MC){
    table = generate_table(N = N, p_a = p_a, p_b = p_b, osn = osn)  # generates table
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

# replace failures in N_LP
for (s in 1:length(scenarios)){
  results_scenarios[[s]][results_scenarios[[s]][,"N_LP"]>10*mean(results_scenarios[[s]][,"N_Chap"]),"N_LP"] = 
    max(results_scenarios[[s]][,"N_Chap"])
}

# organise results for Table 1 and 8
est = numeric()
est_sds = numeric()
est_means = numeric()
for (s in 1:length(scenarios)){
  est_sd = cbind(format(round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1),nsmall=1),round(colSds(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE)/MC^0.5,2))
  est_sd = apply(est_sd,1,paste0, collapse = " (")
  est_sd = paste0(est_sd,")")
  est_sds = cbind(est_sds,est_sd)
  TTESTS = apply(results_scenarios[[s]][,c(2,3,4,8)]-scenarios[[s]]$N,2,t.test)
  est_mean = cbind(format(round(colMeans(results_scenarios[[s]][,c(2,3,4,8)],na.rm=TRUE),1),nsmall=1),as.numeric(unlist(TTESTS)[c(3,14,25,36)]))
  est_mean[as.numeric(est_mean[,2])<=0.001,2] = "***"
  est_mean[as.numeric(est_mean[,2])>0.001&est_mean[,2]<=0.01,2] = "**"
  est_mean[as.numeric(est_mean[,2])>0.01&est_mean[,2]<=0.05,2] = "*"
  est_mean[as.numeric(est_mean[,2])>0.05,2] = ""
  est_mean = apply(est_mean,1,paste0, collapse = "")
  est_means = cbind(est_means,est_mean)
  est = cbind(est,round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1))
}

est_sds = numeric()
for (s in 1:length(scenarios)){
  est_sd  = round(colSds(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1)
  est_sds = cbind(est_sds, est_sd)
}

est_rmses = numeric()
for (s in 1:length(scenarios)){
  est_rmse  = round(colMeans(((results_scenarios[[s]][,c(1,2,3,4,8)]-scenarios[[s]]$N)^2),na.rm=TRUE)^0.5,1)
  est_rmses = cbind(est_rmses, est_rmse)
}

# Results Table 1
t(est_means)
# Results Table 8
cbind(t(est_sds),t(est_rmses))
# end

      
