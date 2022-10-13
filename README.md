# CCFM  
example  
load("TD.rda")  
library(Rcpp)  
library(kernlab)  
sourceCpp("CCM.cpp")  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1)  
optCCMat <- CCM_sel(CCM=CCMat,train_pheno=y.train)  
pheno1 <- predict_pheno(CCMat, optCCMat,"SVM_bag",train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
pheno2 <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=11, chunk=1, method="SVM_bag", CV="ON", num_bootstrap=50)  
