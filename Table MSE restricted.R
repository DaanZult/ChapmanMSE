rm(list=ls())
gc()

set.seed(42)

library(matrixStats)
library(matlib)
library(data.table)
library(brglm2)
library(Rcapture)

# no. iterations
MC = 60000

results_scenarios = list()
scenarios = list()
scenarios[[1]]  = list(Nsources = 3, p_a = 0.5, p_b = 0.4,  p_c = 0.3, N = 100,   osn_AB = 1, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c")
scenarios[[2]]  = list(Nsources = 3, p_a = 0.4, p_b = 0.3,  p_c = 0.2, N = 500,   osn_AB = 1, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c")
scenarios[[3]]  = list(Nsources = 3, p_a = 0.35, p_b = 0.3, p_c = 0.25, N = 10000, osn_AB = 1, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c")

scenarios[[4]]  = list(Nsources = 3, p_a = 0.5, p_b = 0.4,  p_c = 0.3, N = 100,   osn_AB = 1.5, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b")
scenarios[[5]]  = list(Nsources = 3, p_a = 0.4, p_b = 0.3,  p_c = 0.2, N = 500,   osn_AB = 1.5, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b")
scenarios[[6]]  = list(Nsources = 3, p_a = 0.35, p_b = 0.3, p_c = 0.25, N = 10000, osn_AB = 1.5, osn_AC = 1, osn_BC = 1, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b")

scenarios[[7]]  = list(Nsources = 3, p_a = 0.5, p_b = 0.4,  p_c = 0.3, N = 100,   osn_AB = 1.5, osn_AC = 1, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_b:list_c")
scenarios[[8]]  = list(Nsources = 3, p_a = 0.4, p_b = 0.3,  p_c = 0.2, N = 500,   osn_AB = 1.5, osn_AC = 1, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_b:list_c")
scenarios[[9]]  = list(Nsources = 3, p_a = 0.35, p_b = 0.3, p_c = 0.25, N = 10000, osn_AB = 1.5, osn_AC = 1, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_b:list_c")

scenarios[[10]] = list(Nsources = 3, p_a = 0.5, p_b = 0.4,  p_c = 0.3, N = 100,   osn_AB = 1.5, osn_AC = 0.75, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_a:list_c + list_b:list_c")
scenarios[[11]] = list(Nsources = 3, p_a = 0.4, p_b = 0.3,  p_c = 0.2, N = 500,   osn_AB = 1.5, osn_AC = 0.75, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_a:list_c + list_b:list_c")
scenarios[[12]] = list(Nsources = 3, p_a = 0.35, p_b = 0.3, p_c = 0.25, N = 10000, osn_AB = 1.5, osn_AC = 0.75, osn_BC = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_a:list_b + list_a:list_c + list_b:list_c")

scenarios[[13]] = list(Nsources = 4, p_a = 0.4, p_b = 0.35,  p_c = 0.25, p_d = 0.2, N = 20000, osn_AB = 1,   osn_BC = 1,   osn_CD = 1,   osn_AD = 1,   formula = "mysum ~ list_a + list_b + list_c + list_d")
scenarios[[14]] = list(Nsources = 4, p_a = 0.4, p_b = 0.35,  p_c = 0.25, p_d = 0.2, N = 20000, osn_AB = 1.5, osn_AD = 1.5, osn_BC = 0.75, osn_CD = 0.5, formula = "mysum ~ list_a + list_b + list_c + list_d + list_a:list_b + list_b:list_c + list_c:list_d + list_a:list_d + list_a:list_c + list_b:list_d + list_a:list_b:list_c + list_a:list_b:list_d + list_a:list_c:list_d + list_b:list_c:list_d")

generate_table_3sources = function(N = N, Nsources = Nsources, p_a = p_a, p_b = p_b, p_c = p_c, osn_AB = osn_AB, osn_AC = osn_AC, osn_BC = osn_BC){
  ### create dataset from contingency table
  osn_1a = osn_AB
  osn_1b = osn_AB
  
  osn_2a = osn_AC
  osn_2b = osn_AC
  
  osn_3a = osn_BC
  osn_3b = osn_BC
  
  list_a <- factor(c("in","in","in","in","missed","missed","missed","missed"),
                   levels=c("missed","in"))
  
  list_b <- factor(c("in","in","missed","missed","in","in","missed","missed"),
                   levels=c("missed","in"))
  
  list_c <- factor(c("in","missed","in","missed","in","missed","in","missed"),
                   levels=c("missed","in"))
  
  count <- c((p_a*p_b*p_c*N),((p_a)*(p_b)*(1-p_c)*N),
             ((p_a)*(1-p_b)*(p_c)*N),((p_a)*(1-p_b)*(1-p_c)*N),
             ((1-p_a)*(p_b)*(p_c)*N),((1-p_a)*(p_b)*(1-p_c)*N),
             ((1-p_a)*(1-p_b)*(p_c)*N),((1-p_a)*(1-p_b)*(1-p_c)*N))
  
  data_2 <- data.frame(list_a, list_b, list_c, count)
  data_2$offset <- c(log(osn_1b) + log(osn_2b) + log(osn_3b),
                     log(osn_1a), log(osn_2a), log(1),
                     log(osn_3a), log(1), log(1), log(1))
  
  GLM <- glm(formula = count ~ list_a + list_b + list_c,
             family = poisson(link = "log"),
             data = data_2,
             offset = data_2$offset)
  
  summary(GLM)
  GLM$fitted
  data_2$est_countoff<- c(GLM$fitted)
  or_offset<- (GLM$fitted[8]*GLM$fitted[2])/
    (GLM$fitted[4]*GLM$fitted[6])
  or_offset
  
  #BC|A = 0
  ### check to see if offset gave or = osn_3a
  or_offset2<- (GLM$fitted[8]*GLM$fitted[5])/
    (GLM$fitted[6]*GLM$fitted[7])
  or_offset2
  #AC|B = 0
  ### check to see if offset gave or = osn_2a
  or_offset3<- (GLM$fitted[8]*GLM$fitted[3])/
    (GLM$fitted[4]*GLM$fitted[7])
  or_offset3
  #AB|C = 1
  ### check to see if offset gave or = osn_1b
  or_offset<- (GLM$fitted[1]*GLM$fitted[7])/
    (GLM$fitted[3]*GLM$fitted[5])
  or_offset
  #BC|A = 1
  ### check to see if offset gave or = osn_3b
  or_offset2<- (GLM$fitted[1]*GLM$fitted[4])/
    (GLM$fitted[2]*GLM$fitted[3])
  or_offset2
  #AC|B = 1
  ### check to see if offset gave or = osn_2b
  or_offset3<- (GLM$fitted[1]*GLM$fitted[6])/
    (GLM$fitted[5]*GLM$fitted[2])
  or_offset3
  
  ### check the response of the lists again
  ### reponse should be p_a
  list_a_resp <- (GLM$fitted[1]+GLM$fitted[2]+
                    GLM$fitted[3]+GLM$fitted[4])/
    
    (GLM$fitted[1]+GLM$fitted[2]+
       GLM$fitted[3]+GLM$fitted[4]
     +GLM$fitted[5]+GLM$fitted[6]+
       GLM$fitted[7]+GLM$fitted[8])
  list_a_resp
  ### reponse should be p_b
  list_b_resp <- (GLM$fitted[1]+GLM$fitted[2]+
                    GLM$fitted[5]+GLM$fitted[6])/
    
    (GLM$fitted[1]+GLM$fitted[2]+
       GLM$fitted[3]+GLM$fitted[4]+
       GLM$fitted[5]+GLM$fitted[6]+
       GLM$fitted[7]+GLM$fitted[8])
  list_b_resp
  ### reponse should be p_c
  list_c_resp <- (GLM$fitted[1]+GLM$fitted[3]+
                    GLM$fitted[5]+GLM$fitted[7])/
    
    (GLM$fitted[1]+GLM$fitted[2]+
       GLM$fitted[3]+GLM$fitted[4]+
       GLM$fitted[5]+GLM$fitted[6]+
       GLM$fitted[7]+GLM$fitted[8])
  list_c_resp
  
  ### total size
  data_2$total_samp <- N
  data_2_final <- cbind(data_2,data_2$total_samp)
  data_2_final$prob<- data_2_final$est_countoff/data_2_final$total_samp
  ### to store estimates
  total_est_store_2 <- NULL
  total_pop_store <- NULL
  #set.seed(123)
  ### create multinomail samples
  ### create multivariate hypergeometric samples
  data_2_final$mysum <- rmultinom(1, N, data_2_final$est_countoff)
  #data_2_final$mysum <- t(rmvhyper(nn=1, n=count, k=sum(count)))
  data_2_final$mysum <- as.numeric(data_2_final$mysum)
  ### throw away missing cells from data
  data_2<-data_2_final[!(data_2_final$list_a=="missed"&
                           data_2_final$list_b=="missed"&
                           data_2_final$list_c=="missed"),]
  
  table = data_2[c(1:3,5,4,6,7),]
  
  return(table)
}

generate_table_4sources = function(N = N, Nsources = Nsources, p_a = p_a, p_b = p_b, p_c = p_c, p_d = p_d, 
                                   osn_AB=osn_AB, osn_BC=osn_BC, osn_CD=osn_CD, osn_AD=osn_AD){
  osn_1<-osn_AB
  osn_7<-osn_AB
  osn_10<-osn_AB
  osn_19<-osn_AB
  osn_4<-osn_BC
  osn_9<-osn_BC
  osn_22<-osn_BC
  osn_16<-osn_BC
  osn_6<-osn_CD
  osn_24<-osn_CD
  osn_15<-osn_CD
  osn_17<-osn_CD
  osn_3<-osn_AD
  osn_11<-osn_AD
  osn_13<-osn_AD
  osn_21<-osn_AD
  osn_5<-1
  osn_2<-1
  osn_8<-1
  osn_12<-1
  osn_14<-1
  osn_20<-1
  osn_18<-1
  osn_23<-1
  
  list_a <- factor(c("in","in","in","in","in","in","in","in","missed","missed","missed","missed","missed","missed","missed","missed"),
                   levels=c("missed","in"))
  
  list_b <- factor(c("in","in","missed","missed","in","in","missed","missed","in","in","in","in","missed","missed","missed","missed"),
                   levels=c("missed","in"))
  
  list_c <- factor(c("in","in","in","in","missed","missed","missed","missed","in","in","missed","missed","in","in","missed","missed"),
                   levels=c("missed","in"))
  
  list_d <- factor(c("in","missed","in","missed","in","missed","in","missed","in","missed","in","missed","in","missed","in","missed"),
                   levels=c("missed","in"))
  
  count <- c(   p_a  *    p_b  *    p_c  *    p_d  * N,
                (   p_a  *    p_b  *    p_c  * (1-p_d) * N),
                (   p_a  * (1-p_b) *    p_c  *    p_d  * N),
                (   p_a  * (1-p_b) *    p_c  * (1-p_d) * N),
                (   p_a  *    p_b  * (1-p_c) *    p_d  * N),
                (   p_a  *    p_b  * (1-p_c) * (1-p_d) * N),
                (   p_a  * (1-p_b) * (1-p_c) *    p_d  * N),
                (   p_a  * (1-p_b) * (1-p_c) * (1-p_d) * N),
                ((1-p_a) *    p_b  *    p_c  *    p_d  * N),
                ((1-p_a) *    p_b  *    p_c  * (1-p_d) * N),
                ((1-p_a) *    p_b  * (1-p_c) *    p_d  * N),
                ((1-p_a) *    p_b  * (1-p_c) * (1-p_d) * N),
                ((1-p_a) * (1-p_b) *    p_c  *    p_d  * N),
                ((1-p_a) * (1-p_b) *    p_c  * (1-p_d) * N),
                ((1-p_a) * (1-p_b) * (1-p_c) *    p_d  * N),
                ((1-p_a) * (1-p_b) * (1-p_c) * (1-p_d) * N))
  
  data_2 <- data.frame(list_a, list_b, list_c, list_d, count)
  
  data_2$offset <- c(log(osn_1) + log(osn_2) + log(osn_3) +log(osn_4) + log(osn_5) + log(osn_6),
                     log(osn_7) + log(osn_8) + log(osn_9),
                     log(osn_13)+ log(osn_14) + log(osn_15),
                     log(osn_20),log(osn_10)+ log(osn_11) + log(osn_12),
                     log(osn_19), log(osn_21),log(1),log(osn_16) + log(osn_17) + log(osn_18),
                     log(osn_22), log(osn_23),log(1), log(osn_24), log(1), log(1),log(1))
  
  ### model estimated counts from independence model and include offset term
  GLM <- glm(formula = count ~ list_a + list_b + list_c + list_d,
             family = poisson(link = "log"),
             data = data_2,
             offset = data_2$offset)
  
  summary(GLM)
  GLM$fitted
  data_2$est_countoff<- c(GLM$fitted)
  
  data_2$total_samp <- N
  data_2_final <- cbind(data_2,data_2$total_samp)
  data_2_final$prob<- data_2_final$est_countoff/data_2_final$total_samp
  ### to store estimates
  total_est_store_2 <- NULL
  total_pop_store <- NULL
  
  data_2_final$mysum <- rmultinom(1, N, data_2_final$est_countoff)
  data_2_final$mysum <- as.numeric(data_2_final$mysum)
  ### throw away missing cells from data
  missed = data_2_final[(data_2_final$list_a=="missed"&
                           data_2_final$list_b=="missed"&
                           data_2_final$list_c=="missed"&
                           data_2_final$list_d=="missed"),"mysum"]
  data_2<-data_2_final[!(data_2_final$list_a=="missed"&
                           data_2_final$list_b=="missed"&
                           data_2_final$list_c=="missed"&
                           data_2_final$list_d=="missed"),]
  
  table = data_2[c(1,2,5,3,9,6,4,7,10,11,13,8,12,14,15),]
  
  return(table)
}


estimate_MSEs = function(table = table, Nsources = Nsources, formula = formula){
  
  table_save = table
  n=sum(table$mysum)
  
  ### model estimated counts from independence model
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table)
  
  N_ML <- n + exp(summary(GLM)$coefficients[1,1])
  
  # Evans adjustment
  table[,"mysum"] = table[,"mysum"] + 0.5^(Nsources-1)
  
  ### model estimated counts from independence model
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table)
  
  N_Evans <- n + exp(summary(GLM)$coefficients[1,1])
  
  table = table_save
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table, method = "brglmFit",type="MPL_Jeffreys")
  
  N_Firth <- n + exp(summary(GLM)$coefficients[1,1])
  
  ### model estimated counts from independence model
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table, method = "brglmFit", type = "AS_mean")
  
  N_Kosmidis <- n + exp(summary(GLM)$coefficients[1,1])
  
  ### model estimated counts from independence model
  GLM <- try(glm(formula = formula,
                 family = poisson(link = "log"),
                 data = table, method = "brglmFit", type = "correction"), silent = T)
  
  if (class(GLM)[1]=="try-error"){
    N_Cordeiro = N_Firth}
  if (class(GLM)[1]!="try-error"){
    N_Cordeiro <- n + exp(summary(GLM)$coefficients[1,1])
  }
  
  # Rivest adjustment
  if (formula=="mysum ~ list_a + list_b + list_c"){
    table[rowSums(table[,1:Nsources]=="in")==1,"mysum"] = table[rowSums(table[,1:Nsources]=="in")==1,"mysum"] + 1/6
    table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] = table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] + 1/3
  }
  
  if (formula!="mysum ~ list_a + list_b + list_c"&Nsources==3){
    table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] = table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] + 2/3
  }
  
  if (formula!="mysum ~ list_a + list_b + list_c + list_d"&Nsources==4){
    table[rowSums(table[,1:Nsources]=="in")==1,"mysum"] = table[rowSums(table[,1:Nsources]=="in")==1,"mysum"] + 1/8
    table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] = table[rowSums(table[,1:Nsources]=="in")==2,"mysum"] + 1/3
  }
  
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table)
  
  N_Rivest <- n + exp(summary(GLM)$coefficients[1,1])
  
  table = table_save
  
  if (Nsources==3){
    freqs = table$mysum
    X = as.matrix(table[,c(1:Nsources)])
    X[X=="in"]=1
    X[X=="missed"]=0
    X = as.data.frame(X)
    X = t(apply(X,1,as.numeric))
    X = as.data.frame(cbind(1,X))
    colnames(X) = c("mysum", "list_a", "list_b", "list_c")
    formula_X = formula
    formula_X = gsub("list_a:list_b", "I(list_a*list_b)", formula_X)
    formula_X = gsub("list_a:list_c", "I(list_a*list_c)", formula_X)
    formula_X = gsub("list_b:list_c", "I(list_b*list_c)", formula_X)
    X = as.matrix(model.frame(as.formula(formula_X), data = X))
    X=as.matrix(X)
    
    Z=Inverse(t(X)%*%X)%*%t(X)
    z_abc = Z[1,]
    
    z_abc_big0 = z_abc
    z_abc_big0[z_abc_big0>0] = 0
  }
  
  if (Nsources==4){
    freqs = table$mysum
    X = as.matrix(table[,c(1:Nsources)])
    X[X=="in"]=1
    X[X=="missed"]=0
    X = as.data.frame(X)
    X = t(apply(X,1,as.numeric))
    X = as.data.frame(cbind(1,X))
    colnames(X) = c("mysum", "list_a", "list_b", "list_c", "list_d")
    formula_X = formula
    formula_X = gsub(" list_a:list_b ", " I(list_a*list_b) ", formula_X)
    formula_X = gsub(" list_b:list_c ", " I(list_b*list_c) ", formula_X)
    formula_X = gsub(" list_c:list_d ", " I(list_c*list_d) ", formula_X)
    formula_X = gsub(" list_a:list_d ", " I(list_a*list_d) ", formula_X)
    formula_X = gsub(" list_a:list_c ", " I(list_a*list_c) ", formula_X)
    formula_X = gsub(" list_b:list_d ", " I(list_b*list_d) ", formula_X)
    formula_X = gsub("list_a:list_b:list_c", " I(list_a*list_b*list_c) ", formula_X)
    formula_X = gsub("list_a:list_b:list_d", " I(list_a*list_b*list_d) ", formula_X)
    formula_X = gsub("list_a:list_c:list_d", " I(list_a*list_c*list_d) ", formula_X)
    formula_X = gsub("list_b:list_c:list_d", " I(list_b*list_c*list_d) ", formula_X)
    X = as.matrix(model.frame(as.formula(formula_X), data = X))
    X = as.matrix(X)
    
    Z=Inverse(t(X)%*%X)%*%t(X)
    z_abc = Z[1,]
    
    z_abc_big0 = z_abc
    z_abc_big0[z_abc_big0>0] = 0
  }
  
  table$mysum = table$mysum-z_abc_big0
  
  GLM <- glm(formula = formula,
             family = poisson(link = "log"),
             data = table)
  
  N_Chap <- n + exp(summary(GLM)$coefficients[1,1])
  
  return(list(n, N_ML, N_Evans, N_Firth, N_Cordeiro, N_Kosmidis, N_Rivest, N_Chap))
}


for (s in 1:length(scenarios)){
  N = scenarios[[s]]$N
  Nsources = scenarios[[s]]$Nsources
  p_a = scenarios[[s]]$p_a
  p_b = scenarios[[s]]$p_b
  p_c = scenarios[[s]]$p_c
  
  if (Nsources == 4){
    p_d = scenarios[[s]]$p_d
  }
  
  if (Nsources == 3){
    osn_AB <- scenarios[[s]]$osn_AB
    osn_AC <- scenarios[[s]]$osn_AC
    osn_BC <- scenarios[[s]]$osn_BC
  }
  
  if (Nsources == 4){
    p_d = scenarios[[s]]$p_d
    osn_AB <- scenarios[[s]]$osn_AB
    osn_BC <- scenarios[[s]]$osn_BC
    osn_CD <- scenarios[[s]]$osn_CD
    osn_AD <- scenarios[[s]]$osn_AD
  }
  
  formula <- scenarios[[s]]$formula
  
  n = NULL
  N_ML <- NULL
  N_Evans <- NULL
  
  N_Firth <- NULL
  N_Kosmidis <- NULL
  N_Cordeiro <- NULL
  
  N_Rivest <- NULL
  N_Chap <- NULL
  
  for (mc in 1:MC){
    if (Nsources==3){
      table = generate_table_3sources(N = N, Nsources = Nsources, p_a = p_a, p_b = p_b, p_c = p_c, osn_AB = osn_AB, osn_AC = osn_AC, osn_BC = osn_BC)
    }
    
    if (Nsources==4){
      table = generate_table_4sources(N = N, Nsources = Nsources, p_a = p_a, p_b = p_b, p_c = p_c, p_d = p_d, osn_AB = osn_AB, osn_BC = osn_BC, osn_CD = osn_CD, osn_AD = osn_AD)
    }
    
    MSEs  = estimate_MSEs(table, Nsources = Nsources, formula = formula)
    
    n = c(n, MSEs[[1]])
    N_ML       <- c(N_ML,       MSEs[[2]])
    N_Evans    <- c(N_Evans,    MSEs[[3]])
    
    N_Firth    <- c(N_Firth,     MSEs[[4]])
    N_Kosmidis <- c(N_Kosmidis,  MSEs[[5]])
    N_Cordeiro <- c(N_Cordeiro,  MSEs[[6]])
    
    N_Rivest   <- c(N_Rivest,    MSEs[[7]])
    N_Chap     <- c(N_Chap,      MSEs[[8]])
    }
  
  results_scenarios[[s]]=as.matrix(cbind(n, N_ML, N_Evans, 
                                         N_Firth, N_Kosmidis, N_Cordeiro,
                                         N_Rivest, N_Chap))
  print(s)
  print(round(colMeans(results_scenarios[[s]]),1))
}

save.image("C:/Users/daanz/OneDrive/Bureaublad/Werk/Vangst Hervangst/Chapman/Results/Table MSE restricted publicatie.RData")
load("C:/Users/daanz/OneDrive/Bureaublad/Werk/Vangst Hervangst/Chapman/Results/Table MSE restricted publicatie.RData")


for (s in 1:length(scenarios)){
  results_scenarios[[s]][results_scenarios[[s]][,"N_ML"]>10*mean(results_scenarios[[s]][,"N_Chap"]),"N_ML"] = 
    results_scenarios[[s]][results_scenarios[[s]][,"N_ML"]>10*mean(results_scenarios[[s]][,"N_Chap"]),"N_Chap"]
}

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
  est_pval[as.numeric(est_pval[,2])<=0.01,2] = "***"
  est_pval[as.numeric(est_pval[,2])>0.01&est_pval[,2]<=0.05,2] = "**"
  est_pval[as.numeric(est_pval[,2])>0.05&est_pval[,2]<=0.1,2] = "*"
  est_pval[as.numeric(est_pval[,2])>0.1,2] = ""
  est_pval = apply(est_pval,1,paste0, collapse = "")
  est_pvals = cbind(est_pvals,est_pval)
  est = cbind(est,t(t(round(colMeans(results_scenarios[[s]][,c(1,2,3,4,8)],na.rm=TRUE),1))))
}

X = cbind(1:14,
          c(100,500,10000,100,500,10000,100,500,10000,100,500,10000,20000,20000),
          t(est_pvals))
Y = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
Y = gsub(".","&.",Y,fixed=TRUE)
#Y = gsub("0 ","0&. ",Y,fixed=TRUE)
print(Y)

X = cbind(1:14,
          c(100,500,10000,100,500,10000,100,500,10000,100,500,10000,20000,20000),
          t(est_pvals))
MSE_saturated_means = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(MSE_saturated_means)
MSE_saturated_means = gsub(".","&.",MSE_saturated_means,fixed=TRUE)
MSE_saturated_means = gsub("  ","",MSE_saturated_means)
print(MSE_saturated_means)


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
MSE_saturated_sigmas = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(MSE_saturated_sigmas)
MSE_saturated_sigmas = gsub(".","&.",MSE_saturated_sigmas,fixed=TRUE)
print(MSE_saturated_sigmas)

est_rmses = numeric()
for (s in 1:length(scenarios)){
  est_rmse  = format(round(colMeans(((results_scenarios[[s]][,c(1,2,3,4,8)]-scenarios[[s]]$N)^2),na.rm=TRUE)^0.5,1),nsmall=1)
  est_rmses = cbind(est_rmses, est_rmse)
}

est_rmses = est_rmses[-1,]

X = cbind(1:14,
          c(100,500,10000,100,500,10000,100,500,10000,100,500,10000,20000,20000),
          t(est_rmses))
MSE_saturated_rmses = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
print(MSE_saturated_rmses)
MSE_saturated_rmses = gsub(".","&.",MSE_saturated_rmses,fixed=TRUE)
MSE_saturated_rmses = gsub("  ","",MSE_saturated_rmses)
print(MSE_saturated_rmses)

# est_sds_rmses = rbind(est_sds[1,],est_rmses[1,],est_sds[2,],est_rmses[2,],est_sds[3,],est_rmses[3,],est_sds[4,],est_rmses[4,])
# 
# X = cbind(1:7,
#           c(100,100,500,500,10000,10000,100),
#           c(0.5,0.35,0.4,0.25,0.3,0.25,0.15),
#           c(0.2,0.3,0.15,0.2,0.1,0.15,0.15),
#           t(est_sds_rmses))
# MSE_saturated_sigmas_rmses = as.matrix(apply(X, 1, paste0, collapse = " & "),nrow(X)-1,1)
# print(MSE_saturated_sigmas_rmses)
# MSE_saturated_sigmas_rmses = gsub(".","&.",MSE_saturated_sigmas_rmses,fixed=TRUE)
# MSE_saturated_sigmas_rmses = gsub("  ","",MSE_saturated_sigmas_rmses)
# print(MSE_saturated_sigmas_rmses)

MSE_saturated_means
MSE_saturated_sigmas
MSE_saturated_rmses
