
#' @title one_predict_pheno
#' One step finish CCFM
#' @description Users will only get the prediction phenotype of test group.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param train_geno A matrix of the genotype with 0, 1, and 2 for the training group. Here, the column number equals the marker number, and the row number equals the individual number of the train group. Genotype class of the training group can be “big.matrix???.
#' @param test_geno A matrix of the genotype with 0, 1, and 2 for the test group. Here, the column number equals the marker number, which must match the number of markers in the training group, and the row number equals the individual number of the test group. Genotype class of the test group must match the genotype class of the training group.
#' @param CCN The number value is the number of building compressed components.
#' @param chunk If RAM is enough to run, chunk can be 1. But if the RAM is not enough for big data, users can choose an appropriate number to chop the genotype into chunks.
#' @param method The option can be LM, SVM, SVM_bag, LM_Max, SVM_Max, LM_Min, SVM_Min.
#' @param CV If OFF, the argument of method can be LM, SVM, SVM_bag, and If NO and CCN>10 , the argument of method can be LM_Max, SVM_Max, LM_Min, SVM_Min.
#' @param num_bootstrap It's the number of samples in bagging_SVM. The default option is 50.
#' @import RcppEigen kernlab
#' @return Predition_Pheno It is the prediction phenotype of the test group, which order follows the genotype of the test group.
#' @export
# One step get prediction phenotype

one_predict_pheno <- function(train_pheno, train_geno, test_geno, CCN, chunk, method, CV, num_bootstrap=50){
  if ("big.matrix" %in% class(train_geno)){
    CC <- CCmadeBM(bX_train = train_geno@address,bX_test = test_geno@address,by_train = train_pheno, chunk=chunk, cc_num=CCN)
    }else{
      CC <- CCmadeM(y_train = train_pheno, X_train = train_geno, X_test=test_geno, chunk=chunk, cc_num=CCN)
      }
  CCname <- unlist(lapply(1:CCN, function(x){paste0("CC", paste0(x, collapse=""), sep = "")}))
  colnames(CC$CCtrain) <- CCname
  colnames(CC$CCtest) <- CCname
	CCM_train <- CC$CCtrain
	CCM_test <- CC$CCtest
	if(CCN > 10 && CV == "ON"){
	  CV == "ON"
	  }else{
	    CV == "OFF"
	    }
	if(CV == "OFF"){
	  if(method == "LM"){
	    lm = RcppEigen::fastLmPure(X=matrix(cbind(1, CCM_train), ncol= CCN+1), y = train_pheno, method=2)
	    Predition_Pheno <- CCM_test%*%(lm$coefficients[-1])
	    }
	  if(method == "SVM"){
	    svm<-kernlab::ksvm(x = CCM_train, y = train_pheno, kernel="rbfdot")
      Predition_Pheno <- kernlab::predict(svm, CCM_test)
      }
	  if(method == "SVM_bag"){
	    pred <-matrix(0, num_bootstrap, dim(test_geno)[1])
	    for (j in 1:num_bootstrap){
	      train<-cbind(train_pheno, train_geno)
			  train_sample<-train[sample(1:nrow(train), floor(0.2*nrow(train)), replace =T), ]
			  svm_sample<- kernlab::ksvm(x= train_sample[, -1], y= train_sample[, 1], kernel="rbfdot")
			  pred[j, ]<- kernlab::predict(svm_sample, test_geno)
			  }
	    pred_norm <- colMeans(pred, na.rm =T)
	    Predition_Pheno <- (pred_norm-min(pred_norm))/(max(pred_norm)-min(pred_norm))
	    }
	  }else{
	    A.combn  <-  function(n){
	      set <- 0:(2^n-1)
	      comb <- matrix(0, ncol = n, nrow = 2^n)
	      for (i in 1:n){
	        comb[, i] = ifelse((set-rowSums(comb*rep(c(2^((n-1):0)), each=2^n)))/(2^(n-i))>=1, 1, 0)
	        }
	      return(comb)
	      }
	    N<-ncol(CCM_train)
	    K<-ceiling(N/10)
	    foldm <- cut(seq_len(N), breaks=K, labels=FALSE)
      e_best<-vector()
      e_min<-vector()
      for (j in 1:K){
        ind <- which(foldm==j, arr.ind=TRUE)
        ee_best<-vector()
        ee_min<-vector()
        name_set = c(colnames(CCM_train[,ind]))
        n <- length(name_set)
        selectMatrix <- A.combn(n)
        group <- apply(selectMatrix, 1, function(x){name_set[which(x==1)]})
        group <- group[-1]
		    indexes<-list(NULL)
        length(indexes)<-length(group)
        for (m in 1:length(group)) {
          name<-c(unlist(group[m]))
          indexes[m]<-list(name)
			    eee_best<-vector()
			    eee_min<-vector()
			    foldn <- cut(seq_len(dim(CCM_train)[1]), breaks=5, labels=FALSE)
			    for (k in 1:5){
			      index <- which(foldn==k, arr.ind=TRUE)
			     	mod = RcppEigen::fastLmPure(X=as.matrix(cbind(1, CCM_train[-index,name]), ncol=length(name)+1), y= train_pheno[-index], method=2)
			     	if (length(name)==1){
			     	  pred <- CCM_train[index,name]*mod$coefficients[-1]
			     	  }else{
			     	    pred <- CCM_train[index,name]%*%mod$coefficients[-1]
			     	    }
			     	eee_best[k]<- cor(pred, train_pheno[index])
			     	eee_min[k]<- sqrt(mean((pred - train_pheno[index])^2))
			     	}
			    ee_best[m]<-mean(eee_best)
			    ee_min[m]<-mean(eee_min)
			    }
        index_best<- unlist(indexes[which.max(ee_best)])
        index_min<-unlist(indexes[which.min(ee_min)])
        e_best<-c(e_best, index_best)
        e_min<-c(e_min, index_min)
        }
      optCC <- list(best_CC = e_best, min_CC = e_min)
      best_CC <- optCC$best_CC
      min_CC <- optCC$min_CC
      if(method == "LM_Max"){
        lm_best = RcppEigen::fastLmPure(X=cbind(1,CCM_train[,best_CC]), y = train_pheno, method=2)
        Predition_Pheno <- CCM_test[,best_CC]%*%(lm_best$coefficients[-1])
        }
      if(method == "SVM_Max"){
        svm_best<-kernlab::ksvm(x = CCM_train[,best_CC], y = train_pheno, kernel="rbfdot")
        Predition_Pheno <- kernlab::predict(svm_best, CCM_test[,best_CC])
        }
      if(method == "LM_Min"){
        lm_min = RcppEigen::fastLmPure(X=cbind(1,CCM_train[,min_CC]), y = train_pheno, method=2)
        Predition_Pheno <- CCM_test[,min_CC]%*%(lm_min$coefficients[-1])
        }
      if(method == "SVM_Min"){
        svm_min<-kernlab::ksvm(x = CCM_train[,min_CC], y = train_pheno, kernel="rbfdot")
        Predition_Pheno <- kernlab::predict(svm_min, CCM_test[,min_CC])
        }
      }
	return(Predition_Pheno)
	}

#' @title CCM
#' building Compressed component matrix
#' @description building Compressed component matrix for train group and test group.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param train_geno A matrix of the genotype with 0, 1, and 2 for the training group. Here, the column number equals the marker number, and the row number equals the individual number of the train group. Genotype class of the training group can be “big.matrix???.
#' @param test_geno A matrix of the genotype with 0, 1, and 2 for the test group. Here, the column number equals the marker number, which must match the number of markers in the training group, and the row number equals the individual number of the test group. Genotype class of the test group must match the genotype class of the training group.
#' @param CCN The number value is the number of building compressed components.
#' @param chunk If RAM is enough to run, chunk can be 1. But if the RAM is not enough for big data, users can choose an appropriate number to chop the genotype into chunks.
#' @return a list of CCM of train group and test group
#' @export
# building CCM
CCM <- function(train_pheno, train_geno, test_geno, CCN, chunk){
  if ("big.matrix" %in% class(train_geno) ){
    CC <- CCmadeBM(bX_train = train_geno@address,bX_test = test_geno@address,by_train = train_pheno, chunk=chunk, cc_num=CCN)
  }else{
    CC <- CCmadeM(y_train = train_pheno, X_train = train_geno, X_test=test_geno, chunk=chunk, cc_num=CCN)
  }
  CCname <- unlist(lapply(1:CCN, function(x){paste0("CC", paste0(x, collapse=""), sep = "")}))
  colnames(CC$CCtrain) <- CCname
  colnames(CC$CCtest) <- CCname
  CCM_train <- CC$CCtrain
  CCM_test <- CC$CCtest
  return(list(trainCCM=CCM_train,testCCM=CCM_test))
}

#' @title CCM_sel
#' @description Users will select fixed CC from CCM.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param CCM Obtained by CCM function
#' @return A list include CC sets of optimization of prediction ability and CC sets of optimization of prediction error.
#' @export
#'
#'
CCM_sel <- function(CCM,train_pheno){
  CCM_train <- CCM$trainCCM
  CCM_test <- CCM$testCCM
  if(dim(CCM_train)[2] > 10){
    A.combn <- function(n) {
      set <- 0:(2^n - 1)
      comb <- matrix(0, ncol = n, nrow = 2^n)
      for (i in 1:n) {
        comb[, i] = ifelse((set - rowSums(comb * rep(c(2^((n - 1):0)), each = 2^n)))/(2^(n - i)) >= 1, 1, 0)
      }
      return(comb)
    }
    N <- ncol(CCM_train)
    K <- ceiling(N/10)
    foldm <- cut(seq_len(N), breaks = K, labels = FALSE)
    e_best <- vector()
    e_min <- vector()
    for (j in 1:K) {
      ind <- which(foldm == j, arr.ind = TRUE)
      ee_best <- vector()
      ee_min <- vector()
      name_set = c(colnames(CCM_train[, ind]))
      n <- length(name_set)
      selectMatrix <- A.combn(n)
      group <- apply(selectMatrix, 1, function(x) {
        name_set[which(x == 1)]
      })
      group <- group[-1]
      indexes <- list(NULL)
      length(indexes) <- length(group)
      for (m in 1:length(group)) {
        name <- c(unlist(group[m]))
        indexes[m] <- list(name)
        eee_best <- vector()
        eee_min <- vector()
        foldn <- cut(seq_len(dim(CCM_train)[1]), breaks = 5,
                     labels = FALSE)
        for (k in 1:5) {
          index <- which(foldn == k, arr.ind = TRUE)
          mod = RcppEigen::fastLmPure(X = as.matrix(cbind(1,CCM_train[-index, name]), ncol = length(name) +1), 
                                      y = train_pheno[-index], method = 2)
          if (length(name) == 1) {
            pred <- CCM_train[index, name] * mod$coefficients[-1]
          }
          else {
            pred <- CCM_train[index, name] %*% mod$coefficients[-1]
          }
          eee_best[k] <- cor(pred, train_pheno[index])
          eee_min[k] <- sqrt(mean((pred - train_pheno[index])^2))
        }
        ee_best[m] <- mean(eee_best)
        ee_min[m] <- mean(eee_min)
      }
      index_best <- unlist(indexes[which.max(ee_best)])
      index_min <- unlist(indexes[which.min(ee_min)])
      e_best <- c(e_best, index_best)
      e_min <- c(e_min, index_min)
    }
    return(list(max_CC = e_best, min_CC = e_min))
  }
  else{
    return(print("please build CCM of cc_number over 10"))
  }
}

#' @title predict_pheno
#' @description Users will only get the prediction phenotype of test group.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param train_geno A matrix of the genotype with 0, 1, and 2 for the training group. Here, the column number equals the marker number, and the row number equals the individual number of the train group. Genotype class of the training group can be “big.matrix???.
#' @param test_geno A matrix of the genotype with 0, 1, and 2 for the test group. Here, the column number equals the marker number, which must match the number of markers in the training group, and the row number equals the individual number of the test group. Genotype class of the test group must match the genotype class of the training group.
#' @param CCM Obtained by CCM function
#' @param method The option can be LM, SVM, SVM_bag, LM_Max, SVM_Max, LM_Min, SVM_Min.
#' @param num_bootstrap It's the number of samples in bagging_SVM. The default option is 50.
#' @return Predition_Pheno It is the prediction phenotype of the test group, which order follows the genotype of the test group.
#' @export


predict_pheno <- function(method,CCM, train_pheno,test_geno,train_geno,num_bootstrap=50){
  CCM_train <- CCM$trainCCM
  CCM_test <- CCM$testCCM
  if (method == "LM") {
    CCN <- dim(CCM_train)[2]
    lm <- RcppEigen::fastLmPure(X = matrix(cbind(1, CCM_train), ncol = CCN + 1),
                                y = train_pheno,
                                method = 2)
    Predition_Pheno <- CCM_test %*% (lm$coefficients[-1])
  }
  if (method == "SVM") {
    CCN <- dim(CCM_train)[2]
    svm <- kernlab::ksvm(x = CCM_train,
                         y = train_pheno,
                         kernel = "rbfdot")
    Predition_Pheno <- kernlab::predict(svm, CCM_test)
  }
  if (method == "SVM_bag") {
    CCN <- dim(CCM_train)[2]
    pred <- matrix(0, num_bootstrap, dim(test_geno)[1])
    for (j in 1:num_bootstrap) {
      train <- cbind(train_pheno, train_geno)
      train_sample <- train[sample(1:nrow(train), floor(0.2 * nrow(train)), replace = T), ]
      svm_sample <- kernlab::ksvm(x = train_sample[, -1], y = train_sample[, 1], kernel = "rbfdot")
      pred[j, ] <- kernlab::predict(svm_sample, test_geno)
    }
    pred_norm <- colMeans(pred, na.rm = T)
    Predition_Pheno <- (pred_norm - min(pred_norm))/(max(pred_norm) -  min(pred_norm))
  }
  if (method == "LM_Max") {
    opt <- CCM_sel(CCM,train_pheno)
    max_CC <- opt$max_CC
    lm_best <- RcppEigen::fastLmPure(X = cbind(1, CCM_train[, max_CC] ),
                                     y = train_pheno,
                                     method = 2)
    Predition_Pheno <- CCM_test[, max_CC] %*% (lm_best$coefficients[-1])
  }
  if (method == "SVM_Max") {
    opt <- CCM_sel(CCM,train_pheno)
    max_CC <- opt$max_CC
    svm_best <- kernlab::ksvm(x = CCM_train[, max_CC],
                              y = train_pheno,
                              kernel = "rbfdot")
    Predition_Pheno <- kernlab::predict(svm_best, CCM_test[, max_CC])
  }
  if (method == "LM_Min") {
    opt <- CCM_sel(CCM,train_pheno)
    min_CC <- opt$min_CC
    lm_min <- RcppEigen::fastLmPure(X = cbind(1, CCM_train[, min_CC]),
                                    y = train_pheno, method = 2)
    Predition_Pheno <- CCM_test[, min_CC] %*% (lm_min$coefficients[-1])
  }
  if (method == "SVM_Min") {
    opt <- CCM_sel(CCM,train_pheno)
    min_CC <- opt$min_CC
    svm_min <- kernlab::ksvm(x = CCM_train[, min_CC],
                             y = train_pheno, kernel = "rbfdot")
    Predition_Pheno <- kernlab::predict(svm_min, CCM_test[, min_CC])
  }
  return(Predition_Pheno)
}

#' @title pre_test
#' One step finish CCFM
#' @description Users will only get the prediction phenotype of test group.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param train_geno A matrix of the genotype with 0, 1, and 2 for the training group. Here, the column number equals the marker number, and the row number equals the individual number of the train group. Genotype class of the training group can be “big.matrix???.
#' @param pre_test_num It's the number of samples in pre-experiment.
#' @return The number of times the method with the highest prediction accuracy appeared in the pre-experiment.
#' @export


pre_test <- function(train_pheno,train_geno,pre_test_num){
  accuracy <- matrix(nrow = pre_test_num, ncol = 11)
  colnames(accuracy) <- c("test_num","LM5","LM6","SVM5","SVM6","Best_LM",
                           "Min_LM","Best_SVM","Min_SVM","svmbag10","svmbag20")
  method <- vector(length = pre_test_num)
  for(hoop in 1:pre_test_num){
    Num <- nrow(train_geno)
    Indexes <- sample(1:Num,Num*0.8,replace=F)
    y_train_data <- train_pheno[Indexes]
    X_train_data <- train_geno[Indexes,]
    y_test_data <- train_pheno[-Indexes]
    X_test_data <- train_geno[-Indexes,]
    CC <- CCmadeM(y_train = y_train_data, X_train = X_train_data, X_test=X_test_data, chunk=1, cc_num=20)
    CCname <- unlist(lapply(1:20, function(x){paste0("CC", paste0(x, collapse=""), sep = "")}))
    colnames(CC$CCtrain) <- CCname
    colnames(CC$CCtest) <- CCname
    XX1 <- CC$CCtrain
    XX1.1 <- CC$CCtest
    CCR_LM <- function(n,y.train, X.train, X.test, chunk=NULL, cc_num=10){
      lm1 = RcppEigen::fastLmPure(X=matrix(cbind(1, XX1[, 1:n]), ncol=n+1), y=y.train, method=2)
      Predition_Pheno <- XX1.1[, 1:n]%*%(lm1$coefficients[-1])
      return(Predition_Pheno)
    }
    CCR_SVM <- function(n,y.train, X.train, X.test, chunk=NULL, cc_num=10){
      mm<-kernlab::ksvm(x= XX1[,1:n], y=y.train, kernel="rbfdot")
      Predition_Pheno <- kernlab::predict(mm, XX1.1[,1:n])
      return(Predition_Pheno)
    }
    svmbagging <- function (y.train=y.train, X.train=XX1, X.test=XX1.1, num_bootstrap=50) {
      normalize<-function(x){return((x-min(x))/(max(x)-min(x)))}
      haha<-matrix(0, num_bootstrap, dim(X.test)[1])
      for (j in 1:num_bootstrap){
        train1<-cbind(y.train, X.train)
        train2<-train1[sample(1:nrow(train1), floor(0.2*nrow(train1)), replace =T), ]
        mm<- kernlab::ksvm(x= train2[, -1], y= train2[, 1], kernel="rbfdot")
        haha[j, ]<- kernlab::predict(mm, X.test)
      }
      return(normalize(colMeans(haha, na.rm =T)))
    }
    chunked_optimization <- function(y.train = y.train, XX1 = XX1){
      A.combn  <-  function(n){
        set <- 0:(2^n-1)
        comb <- matrix(0, ncol = n, nrow = 2^n)
        for (i in 1:n){
          comb[, i] = ifelse((set-rowSums(comb*rep(c(2^((n-1):0)), each=2^n)))/(2^(n-i))>=1, 1, 0)
        }
        return(comb)
      }
      N<-ncol(XX1)
      KK<-ceiling(N/10)
      foldm <- cut(seq_len(N), breaks=KK, labels=FALSE)
      e1<-vector();e2<-vector()
      for (j in 1:KK){
        ind <- which(foldm==j, arr.ind=TRUE)
        ee1<-vector();ee2<-vector()
        c = c(colnames(XX1[,ind]))
        selectMatrix <- A.combn(length(c))
        hehe<-apply(selectMatrix, 1, function(x){c[which(x==1)]})
        hehe<-hehe[-1]
        indexes<-list(NULL)
        length(indexes)<-length(hehe)
        for (m in 1:length(hehe)) {
          name<-c(unlist(hehe[m]))
          indexes[m]<-list(name)
          eee1<-vector();eee2<-vector()
          foldn <- cut(seq_len(dim(XX1)[1]), breaks=5, labels=FALSE)
          for (k in 1:5){
            ind1 <- which(foldn==k, arr.ind=TRUE)
            mod = RcppEigen::fastLmPure(X=as.matrix(cbind(1, XX1[-ind1,name]), ncol=length(name)+1), y= y.train[-ind1], method=2)
            if (length(name)==1){pred <- XX1[ind1,name]*mod$coefficients[-1]}else{
              pred <- XX1[ind1,name]%*%mod$coefficients[-1]}
            eee1[k]<- cor(pred, y.train[ind1])
            eee2[k]<- sqrt(mean((pred - y.train[ind1])^2))
          }
          ee1[m]<-mean(eee1)
          ee2[m]<-mean(eee2)
        }
        index1<- unlist(indexes[which.max(ee1)])
        index2<-unlist(indexes[which.min(ee2)])
        e1<-c(e1, index1)
        e2<-c(e2, index2)
      }
      return(list(best_CC = e1, min_CC = e2))
    }

    optCC <-chunked_optimization(y.train = y_train_data,XX1 = XX1)
    index <- optCC$best_CC
    index1 <- optCC$min_CC

    Best_LM <- RcppEigen::fastLmPure(X=cbind(1, XX1[, index]), y=y_train_data, method=2)
    Min_LM <- RcppEigen::fastLmPure(X=cbind(1, XX1[, index1]), y=y_train_data, method=2)
    Best_SVM <- kernlab::ksvm(x= XX1[, index], y=y_train_data, kernel="rbfdot")
    Min_SVM <-  kernlab::ksvm(x= XX1[, index1], y=y_train_data, kernel="rbfdot")
    CCR_LM5 <- CCR_LM(5,y_train_data, X_train_data, X_test_data, chunk=NULL, cc_num=10)
    CCR_LM6 <- CCR_LM(6,y_train_data, X_train_data, X_test_data, chunk=NULL, cc_num=10)
    CCR_SVM5 <- CCR_SVM(5,y_train_data, X_train_data, X_test_data, chunk=NULL, cc_num=10)
    CCR_SVM6 <- CCR_SVM(6,y_train_data, X_train_data, X_test_data, chunk=NULL, cc_num=10)
    CCR_Best_LM <- XX1.1[, index]%*%(Best_LM$coefficients[-1])
    CCR_Min_LM <- XX1.1[, index1]%*%(Min_LM$coefficients[-1])
    CCR_Best_SVM <- kernlab::predict(Best_SVM, XX1.1[, index])
    CCR_Min_SVM <- kernlab::predict(Min_SVM, XX1.1[, index1])
    svmbag10 <-svmbagging(y.train=y_train_data, X.train=XX1[,1:10], X.test=XX1.1[,1:10], num_bootstrap=50)
    svmbag20 <-svmbagging(y.train=y_train_data, X.train=XX1[,1:20], X.test=XX1.1[,1:20], num_bootstrap=50)
    Ind <- length(y_test_data)
    test.pred<- matrix(NA,nrow = Ind + 1,ncol =11)
    test.pred[1:Ind,1] <- y_test_data
    test.pred[1:Ind,2] <- CCR_LM5
    test.pred[1:Ind,3] <- CCR_LM6
    test.pred[1:Ind,4] <- CCR_SVM5
    test.pred[1:Ind,5] <- CCR_SVM6
    test.pred[1:Ind,6] <- CCR_Best_LM
    test.pred[1:Ind,7] <- CCR_Min_LM
    test.pred[1:Ind,8] <- CCR_Best_SVM
    test.pred[1:Ind,9] <- CCR_Min_SVM
    test.pred[1:Ind,10] <- svmbag10
    test.pred[1:Ind,11] <- svmbag20
    colnames(test.pred) <- c("y","LM5","LM6","SVM5","SVM6","Best_LM",
                             "Min_LM","Best_SVM","Min_SVM","svmbag10","svmbag20")
    for (i in 1:11) {
      test.pred[Ind+1,i] <- cor(y_test_data,test.pred[1:Ind,i])
    }
    accuracy[hoop,] <- test.pred[Ind+1,]
    accuracy[hoop,1] <- hoop
    method[hoop] <- names(which.max(accuracy[hoop,-1]))
  }
  method_time <- table(method)
  return(method_time)
  }

