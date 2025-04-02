
# CCFM
Compressed Components with Flexible Modeling for genomic prediction
## example  
### prepared the data
library(CCFM)  
library(bigmemory)  
fpath <- system.file("extdata", "ID1000.RDATA", package="CCFM", mustWork = TRUE)  
load(fpath)  
ind <- dim(Xtt)[1]  
set.seed(100)  
testIndexes <- sample(1:ind,size = ind*0.2,replace = FALSE)  
y.test <- Xtt[testIndexes,1]  
X.test <- Xtt[testIndexes,-1]  
y.train <- Xtt[-testIndexes,1]  
X.train <- Xtt[-testIndexes,-1]  
BX.train <- as.big.matrix(X.train)  
BX.test <- as.big.matrix(X.test)

### biuld K=20 CCM  
#### the type of genotype is 'matrix'
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1)  
#### the type of genotype is 'big.matrix'
CCMat <- CCM(train_pheno=y.train, train_geno=BX.train, test_geno=BX.test, CCN=20, chunk=1) 

### LM5 predict phenotype  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=5, chunk=1)  
pheno <- predict_pheno(method='LM',CCM=CCMat, train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
### OR  
phenotype <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=5, chunk=1, method="LM", CV="OFF")  

### SVM5 predict phenotype  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=5, chunk=1)  
pheno <- predict_pheno(method='SVM',CCM=CCMat, train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
### OR   
phenotype <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=5, chunk=1, method="SVM", CV="OFF")  

### SVM_bag10 predict phenotype  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=10, chunk=1)  
pheno <- predict_pheno(method='SVM_bag',CCM=CCMat, train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
### OR   
phenotype <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=5, chunk=1, method="SVM_bag", CV="OFF")  
#### there are some difference between pheno and phenotype because of sampling.  

### LM_Max predict phenotype  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1)  
pheno <- predict_pheno(method='LM_Max',CCM=CCMat, train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
### OR   
phenotype <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1, method="LM_Max", CV="ON")  

### SVM_Min predict phenotype  
CCMat <- CCM(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1)  
pheno <- predict_pheno(method='SVM_Min',CCM=CCMat, train_pheno=y.train, train_geno=X.train, test_geno=X.test)  
### OR   
phenotype <- one_predict_pheno(train_pheno=y.train, train_geno=X.train, test_geno=X.test, CCN=20, chunk=1, method="SVM_Min", CV="ON")  

### selecting best method for test data  
name <- pre_test(train_pheno = y.train ,train_geno =X.train ,pre_test_num = 10)  
The method with the highest numerical value is the optimal method




