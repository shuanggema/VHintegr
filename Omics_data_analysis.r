##Part 1 Vertical integrative analysis of multi-omics data
  load("vertical integrative analysis of multi-omics data.RData")

  library(PMA)
  library(boot)
  
  #define function
  vertical_integrative_analysis<-function(GE,CNV,y)
     {
        #spca
        cv_GE        <- SPC.cv(GE,trace = FALSE)                                             # cv choose parameters
        spc_GE       <- SPC(GE, sumabsv=cv_GE$bestsumabsv1se, K=10,trace=FALSE)              # conduct sparse PCA for GE get K=10 pc
        pc_GE        <- spc_GE$u                                                             # data GE pcs
  
        cv_CNV       <- SPC.cv(CNV,trace = FALSE)                                 
        spc_CNV      <- SPC(CNV, sumabsv=cv_CNV$bestsumabsv1se, K=10,trace=FALSE)   
        pc_CNV        <- spc_CNV$u
  
        GE_CNV       <- lm(CNV~scale(pc_GE))$residuals 
        cv_GECNV     <- SPC.cv(GE_CNV,trace = FALSE)
        spc_GECNV    <- SPC(GE_CNV, sumabsv=cv_GECNV$bestsumabsv1se, K=10,trace=FALSE)
        pc_GECNV     <- spc_GECNV$u
  
        #logistic regression
        x<-cbind(pc_GE,pc_GECNV)  #for model A1 
        #x<-cbind(pc_GE,pc_CNV)   #for model A2
        #x<-pc_GE                 #for model A3
        #x<-pc_CNV                #for model A4
        result_mean         <- vector()
  
        for(loop in 1:200)
            {
              #set.seed(loop)
              n<-nrow(y);train_pos<-sample(1:n,3/4*n);test_pos<-setdiff(1:n,train_pos)         # split the data to training set and test set
              B=data.frame(X=x,y)
              fit=glm(y~.,data=B,family="binomial",subset=train_pos,control=list(maxit=200))
              pred=predict(fit,B[test_pos,],type="response")
              result=rep(0,length(test_pos));result[which(pred>.5)]<-1                         
              result_mean<-c(result_mean,mean(result==y[test_pos,]))                           # calculate the corrected prediction ratio (CPR)
            }
  
        return(result_mean=result_mean)
    }

  # vertical analysis of bipolar disorder
  CPR_bipolar<-vertical_integrative_analysis(E_BN,C_BN,y_BN)                                   #the corrected prediction ratios (CPR) of bipolar disease
  mean(CPR_bipolar)                                                                            #the mean of corrected prediction ratios (CPR) of bipolar disease

  # vertical analysis of bipolar disorder
  CPR_schizophrenia<-vertical_integrative_analysis(E_SN,C_SN,y_SN)                             #the corrected prediction ratios (CPR) of schizophrenia
  mean(CPR_schizophrenia)                                                                      #the mean of corrected prediction ratios (CPR) of schizophrenia







##Part 2 Horizontal integrative analysis for disease marker identification
#Take GE for example
   load("Horizontal integrative analysis for disease marker identification.RData")
   source("IntegrLogit.r")

   xi = 0                                             # for B1 model (magnitude-based shrinkage penalty)
   xi=0.0001                                          # for B2 model (sign-based shrinkage penalty)

   select_ge <- matrix(0, nrow = 3, ncol=3)
   colnames(select_ge) <- c("B", "S", "Overlap")
   for (i in 1:3)
    {
        set.seed(20142015)
        x <- list(X[[i]], X[[3+i]])
        nlambda1 =30
        lambda_ratio = 0.05
        lambda2_seq = c(0, 0.01, 0.1, 0.5, 0.7,1)
        S = cbind(sign(cor(x[[1]], y[[1]])), sign(cor(x[[2]], y[[2]])))
        FIT <- IntegrLogit(x, y, S, nlambda1, lambda_ratio, lambda2_seq, xi)
        select_ge[i, 1:2] <- colSums(FIT$beta_cv[-1, ] !=0)
        select_ge[i, 3] <- sum((FIT$beta_cv[-1, 1]!=0) * (FIT$beta_cv[-1, 2]!=0))
    }
   
   select_ge
   
   
  
   
    
   
## Part 3 Horizontal integrative analysis of gene expression-CNV regulations
   load("Horizontal integrative analysis of gene expression-CNV regulations.RData")
   source("IntegrLm.r")
   
   nlambda1 =30
   lambda_ratio = 0.05
   xi = 0                                             # for C1 model (magnitude-based shrinkage penalty)
   xi=0.0000001                                          # for C2 model (sign-based shrinkage penalty)
   
   # for pathway i
   i = 1  #i=2; i=3
   eta_Bipolar <- matrix(0, nrow = ncol(bipolar_cnv[[i]]), ncol = ncol(bipolar_ge[[i]]))
   eta_Schizophrenia <- matrix(0, nrow = ncol(schizophrenia_cnv[[i]]), ncol = ncol(schizophrenia_ge[[i]]))
   
   x = list(bipolar_cnv[[i]], schizophrenia_cnv[[i]])
   for (j in 1:ncol(bipolar_ge[[i]]))	
     {
       y = list(bipolar_ge[[i]][, j],schizophrenia_ge[[i]][, j]) 
       lambda2_seq = c(0,0.01, 0.1, 0.5, 1)
       S = cbind(sign(cor(x[[1]], y[[1]])), sign(cor(x[[2]], y[[2]])))
       FIT <- IntegrLm(x, y, S, nlambda1, lambda_ratio, lambda2_seq, xi)
       eta_Bipolar[, j] = FIT$beta_cv[-1, 1]        
       eta_Schizophrenia[, j] = FIT$beta_cv[-1, 2]
     }
  
   eta_Bipolar                                         #coefficient matrix of bipolar disorder
   eta_Schizophrenia                                   #coefficient matrix of schizophrenia
   
   
   
   
   