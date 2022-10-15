
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
#' @param num_bootstrap It's the number of samples in bagging. The default option is 50.
#' @import RcppEigen kernlab
#' @return Predition_Pheno It is the prediction phenotype of the test group, which order follows the genotype of the test group.
#' @export

# One step get prediction phenotype
one_predict_pheno <- function(train_pheno, train_geno, test_geno, CCN, chunk, method, CV, num_bootstrap=50){
	set.seed(100)
  if (class(train_geno)[1] == "big.matrix"){
	CC <- CCmadeBM(bX_train = train_geno@address,bX_test = test_geno@address,by_train = train_pheno, chunk=chunk, cc_num=CCN)
	}else{
	CC <- CCmadeM(y_train = train_pheno, X_train = train_geno, X_test=test_geno, chunk=chunk, cc_num=CCN)}
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
        Predition_Pheno <- predict(svm, CCM_test)
		}
		if(method == "SVM_bag"){
		pred <-matrix(0, num_bootstrap, dim(test_geno)[1])
		for (j in 1:num_bootstrap){
			train<-cbind(train_pheno, train_geno)
			train_sample<-train[sample(1:nrow(train), floor(0.2*nrow(train)), replace =T), ]
			svm_sample<- kernlab::ksvm(x= train_sample[, -1], y= train_sample[, 1], kernel="rbfdot")
			pred[j, ]<- predict(svm_sample, test_geno)
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
        Predition_Pheno <- predict(svm_best, CCM_test[,best_CC])
		}
	if(method == "LM_Min"){
		lm_min = RcppEigen::fastLmPure(X=cbind(1,CCM_train[,min_CC]), y = train_pheno, method=2)
        Predition_Pheno <- CCM_test[,min_CC]%*%(lm_min$coefficients[-1])
		}
	if(method == "SVM_Min"){
	    svm_min<-kernlab::ksvm(x = CCM_train[,min_CC], y = train_pheno, kernel="rbfdot")
        Predition_Pheno <- predict(svm_min, CCM_test[,min_CC])
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
#'

# building CCM
CCM <- function(train_pheno, train_geno, test_geno, CCN, chunk){
   if (class(train_geno)[1] == "big.matrix"){
	CC <- CCmadeBM(bX_train = train_geno@address,bX_test = test_geno@address,by_train = train_pheno, chunk=chunk, cc_num=CCN)
	}else{
	CC <- CCmadeM(y_train = train_pheno, X_train = train_geno, X_test=test_geno, chunk=chunk, cc_num=CCN)}
  CCname <- unlist(lapply(1:CCN, function(x){paste0("CC", paste0(x, collapse=""), sep = "")}))
  colnames(CC$CCtrain) <- CCname
  colnames(CC$CCtest) <- CCname
  CCM_train <- CC$CCtrain
  CCM_test <- CC$CCtest
  return(list(trainCCM=CCM_train,testCCM=CCM_test))
}

#' @title CCM_sel
#' selecting appropriate CC from CCM
#' @description selecting appropriate CC for building new CCM from original CCM by cross validation.
#'
#' @param CCM By function CCM building CCM of the genotype of train and test group.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @return a list of CCM name index contain CC combination of maximum prediction accuracy and minimum error
#' @export
#'

# selecting approximate CC number
CCM_sel <- function(CCM,train_pheno){
  set.seed(100)
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
          mod = RcppEigen::fastLmPure(X = as.matrix(cbind(1,
                                                          CCM_train[-index, name]), ncol = length(name) +
                                                      1), y = train_pheno[-index], method = 2)
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

#' @title CCM_sel
#' selecting appropriate CC from CCM
#' @description selecting appropriate CC for building new CCM from original CCM by cross validation.
#' @param train_pheno A column is the phenotype of the training group, and the order of phenotype follows the genotype of the training group.
#' @param train_geno A matrix of the genotype with 0, 1, and 2 for the training group. Here, the column number equals the marker number, and the row number equals the individual number of the train group. Genotype class of the training group can be “big.matrix???.
#' @param test_geno A matrix of the genotype with 0, 1, and 2 for the test group. Here, the column number equals the marker number, which must match the number of markers in the training group, and the row number equals the individual number of the test group. Genotype class of the test group must match the genotype class of the training group.
#' @param CCM By function CCM building CCM of the genotype of train and test group.
#' @param opt a list of CCM name index contain CC combination of maximum prediction accuracy and minimum error.
#' @param method The option can be LM, SVM, SVM_bag, LM_Max, SVM_Max, LM_Min, SVM_Min.
#' @param num_bootstrap It's the number of samples in bagging. The default option is 50.
#' @return the prediction phenotype of the test group, which order follows the genotype of the test group.
#' @export
#'

#last step
predict_pheno <- function(CCM, opt, method,train_pheno,test_geno,train_geno,num_bootstrap=50){
  set.seed(100)
  CCM_train <- CCM$trainCCM
  CCM_test <- CCM$testCCM
  CCN <- dim(CCM_train)[2]
  if (method == "LM") {
    lm <- RcppEigen::fastLmPure(X = matrix(cbind(1, CCM_train), ncol = CCN + 1),
                                y = train_pheno,
                                method = 2)
    Predition_Pheno <- CCM_test %*% (lm$coefficients[-1])
  }
  if (method == "SVM") {
    svm <- kernlab::ksvm(x = CCM_train,
                         y = train_pheno,
                         kernel = "rbfdot")
    Predition_Pheno <- predict(svm, CCM_test)
  }
  if (method == "SVM_bag") {
    pred <- matrix(0, num_bootstrap, dim(test_geno)[1])
    for (j in 1:num_bootstrap) {
      train <- cbind(train_pheno, train_geno)
      train_sample <- train[sample(1:nrow(train), floor(0.2 * nrow(train)), replace = T), ]
      svm_sample <- kernlab::ksvm(x = train_sample[, -1], y = train_sample[, 1], kernel = "rbfdot")
      pred[j, ] <- predict(svm_sample, test_geno)
    }
    pred_norm <- colMeans(pred, na.rm = T)
    Predition_Pheno <- (pred_norm - min(pred_norm))/(max(pred_norm) -  min(pred_norm))
  }
  max_CC <- opt$max_CC
  min_CC <- opt$min_CC
  if (method == "LM_Max") {
    lm_best <- RcppEigen::fastLmPure(X = cbind(1, CCM_train[, max_CC]),
                                     y = train_pheno,
                                     method = 2)
    Predition_Pheno <- CCM_test[, max_CC] %*% (lm_best$coefficients[-1])
  }
  if (method == "SVM_Max") {
    svm_best <- kernlab::ksvm(x = CCM_train[, max_CC],
                              y = train_pheno,
                              kernel = "rbfdot")
    Predition_Pheno <- predict(svm_best, CCM_test[, max_CC])
  }
  if (method == "LM_Min") {
    lm_min <- RcppEigen::fastLmPure(X = cbind(1, CCM_train[, min_CC]),
                                    y = train_pheno, method = 2)
    Predition_Pheno <- CCM_test[, min_CC] %*% (lm_min$coefficients[-1])
  }
  if (method == "SVM_Min") {
    svm_min <- kernlab::ksvm(x = CCM_train[, min_CC],
                             y = train_pheno, kernel = "rbfdot")
    Predition_Pheno <- predict(svm_min, CCM_test[, min_CC])
  }
  return(Predition_Pheno)
}



