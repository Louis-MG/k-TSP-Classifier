library(switchBox)

############################################
#
#   TUTORIAL
#
############################################

### LOADING DATA

### Load gene expression data for the test set
data(trainingData)

### Show the class of the ``matTesting'' object
class(matTraining)

### Show the dimentions of the ``matTesting'' matrix
dim(matTraining)

### Show the first  10 sample names of the ``matTest'' matrix
head(colnames(matTraining), n=10)


### ENTRAINEMENT

### Train a classifier using default filtering function based on the Wilcoxon test
classifier <- SWAP.KTSP.Train(matTraining, trainingGroup)

#scores = SWAP.Calculate.BasicTSPScores(testingGroup, matTesting)

classifier

### Apply the classifier to the TRAINING set using default decision rule
trainingPrediction <- SWAP.KTSP.Classify(matTraining, classifier)


### Resubstitution performance in the TRAINING set
### Define a "positive" test result if needed
table(trainingPrediction, trainingGroup)


table(trainingPrediction, trainingGroup)

############################################
#
#   NOS DONNEES
#
############################################

### PREPARE DATA

#setwd("C:/Users/lmgue/OneDrive/M2/PROJET/")
setwd("/home/lmgueguen/Documents/M2/PROJET/")

df = read.csv("projet_data/minusRef_fc2.csv", header = TRUE, dec = ",")
dim(df)
df = unique.data.frame(df)
dim(df)
df = df[,2:58] #removes first column

colnames(df)[1] <- "Genes"
unknown = which(df['Genes'] == '?') #removes genes with "?" as names 
df = df[-unknown,]

df$Genes <- make.names(df$Genes, unique = TRUE)
rownames(df) <- df$Genes
df <- subset(df, select = -c(Genes))

col_names = colnames(df)
col_names = gsub('^SLS.*', 'nonallergique', col_names) #uses regex to change lines starting with SLS (followed by any character any number of times), by irritant
col_names = gsub('^[^[:lower:]].*' ,'allergique' , col_names) #changes lines starting NOT by a lower character and followed by any any number of times by allergique
colnames(df) <- col_names
colnames(df)
dim(df)

### MATRICE

# prepare train and test matrix, 80% and 20% of the data set (respectively 46 and 11 samples)
trainmat = as.matrix(df[,1:49])
rownames(trainmat) <- rownames(df)
testmat = as.matrix(df[,50:56])
rownames(testmat) <- rownames(df)

### ENTRAINEMENT

classifier <- SWAP.KTSP.Train(trainmat, as.factor(colnames(df)[1:49]), krange = c(7))

classifier

# $name
# [1] "7TSPs"
# 
# $TSPs
# [,1]        [,2]      
# [1,] "CLDN23"    "GBP1"    
# [2,] "PLA2G4F"   "PARP9.2" 
# [3,] "SCCPDH"    "PLSCR1.1"
# [4,] "GAL3ST4.1" "CDKN2D"  
# [5,] "LY6K"      "NMI"     
# [6,] "ANKK1"     "DTX3L.1" 
# [7,] "GSTA3"     "ERAP2.1" 
# 
# $score
# [1] 1.0000076 1.0000058 1.0000050 1.0000048 1.0000044 1.0000035 0.9743641
# 
# $labels
# [1] "allergique"     "non-allergique"


prediction1 = SWAP.KTSP.Classify(trainmat, classifier)

reality1 <- as.factor(colnames(df)[1:49])

table(prediction1, reality1)

prediction2 = SWAP.KTSP.Classify(testmat, classifier)

reality2 <- as.factor(colnames(df)[50:56])

table(prediction2, reality2)

                        # reality2
# prediction2      allergique non-allergique
# allergique              6              0
# non-allergique          0              1

library(pROC)
auc(reality2, order(prediction2, decreasing = TRUE)) 
#when setting the decreasing = TRUE, auc of 1 and if the parameter is set on default (decreasing = FALSE), auc = 0.667
#this should be looked at before further exploration


############################################
#
#   CROSS-VALIDATION 
#
############################################


### KFOLD STRATIFIE

#la premi?re fonciton s?pare le jeu de donnles en dataframes par cat?gories, en calculant la proportion du jeu de donn?es initiale pour chaque cat?gorie

sep_by_class = function(df, axis = 2) {
  ###   df a dataframe 
  ###   axis is 1 for samples as columns or 0 for samples as rows
  list_prop <- list() #list of the proportion of each class
  if (axis == 2) {
    class_labels <- unique(colnames(df))
    for (i in 1:(length(class_labels))) {
      prop_class <- length(which(colnames(df) == class_labels[i])) / dim(df)[2]
      list_prop <- append(list_prop, prop_class)
    }
    S3_object <- list(class_labels, list_prop) #initialise S3 object with the proportions for each class
    for (i in 1:length(class_labels)) {
      S3_object <- append(S3_object, list(subset(df, select = which(colnames(df) == class_labels[i])) )) #adds a separated dataframe for each class
    }
    names(S3_object) <- as.list(c("labels", "proportion", class_labels)) #renames the S3 object attributes
  } else {
    class_labels <- unique(rownames(df))
    print(class_labels)
    for (i in 1:length(class_labels)) {
      prop_class <- length(which(rownames(df) == class_labels[i])) / dim(df)[1]
      list_prop <- append(list_prop, prop_class)
    }
    print(list_prop)
    S3_object <- list(class_labels, list_prop) #initialise S3 object with the proportions for each class
    for (i in 1:length(class_labels)) {
      S3_object <- append(S3_object, list(subset(df, select = which(rownames(df) == class_labels[i])) )) #adds a separated dataframe for each class
    }
    names(S3_object) <- as.list(c("labels", "proportion", class_labels)) #renames the S3 object attributes
  }
  for (i in S3_object$labels) {
   colnames(S3_object[[i]]) <- rep(i, dim(S3_object[[i]])[2])
  }
  return(S3_object) #returns S3 class object with proportions for each class, dataframes for the data from each class, in the same order 
}

#does not woooooooooooooork
rename_ <- function(x, i) {
  col_names <- colnames(x) 
  col_names <- gsub(paste('^X*', i[1], sep = ""), i[1], col_names)
  col_names <- gsub(paste('^X*', i[2], sep = ""), i[2], col_names)
  return(col_names)
}

### add if else conditions for rows and columns with arg axis

strat_kfold = function(df, k, axis = 2) {
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  if (axis == 2) {
    data_by_class <- sep_by_class(df, axis = axis) #splits dataframe by class in sub dataframe as S3 class object attributes
    sets_index = list() #list of indexes for the k subsets for kfold cv
    sets_index <- lapply(data_by_class[3:length(data_by_class)], function(x) split(sample(c(1:dim(x)[2]), dim(x)[2]), sort(c(1:dim(x)[2])%%k +1))) #list, for each class there are k lists containing the indexes 
    names(sets_index) <- as.list(data_by_class$labels)
    data_by_class_by_fold <- lapply(data_by_class$labels, function(x) lapply(sets_index[[x]], function(y) data_by_class[[x]][,y]))#here the indexes are used to make a list of dataframes each correspondinfg to a fold
    names(data_by_class_by_fold) <- as.list(data_by_class$labels)
    folds <- list()
    folds_names <- list()
    for (i in c(1:k)) {
      folds <- append(folds, list(as.data.frame(lapply(data_by_class$labels, function(x) data_by_class_by_fold[[x]][i]) )) ) #adds class1 and class2 datasets together to for k sub datasets with equal proportions of each class
      colnames(folds) <- col_names
      folds_names <- append(folds_names, paste("fold", i, sep = "")) 
    }
    names(folds) <- folds_names
  } #else {  #### the next lines do not work
  #   data_by_class <- sep_by_class(df, axis = axis) #splits dataframe by class in sub dataframe as S3 class object attributes
  #   sets_index = list() #list of indexes for the k subsets for kfold cv
  #   sets_index <- lapply(data_by_class[3:length(data_by_class)], function(x) split(sample(c(1:dim(x)[1]), dim(x)[1]), sort(c(1:dim(x)[1])%%k +1))) #list, for each class there are k lists containing the indexes 
  #   names(sets_index) <- as.list(data_by_class$labels)
  #   data_by_class_by_fold <- lapply(data_by_class$labels, function(x) lapply(sets_index[[x]], function(y) data_by_class[[x]][y,]))#here the indexes are used to make a list of dataframes each correspondinfg to a fold
  #   names(data_by_class_by_fold) <- as.list(data_by_class$labels) #pas besoin visiblement
  #   folds <- list()
  #   folds_names <- list()
  #   for (i in c(1:k)) {
  #     folds <- append(folds, list(as.data.frame(lapply(data_by_class$labels, function(x) data_by_class_by_fold[[x]][i]) )) )
  #     folds_names <- append(folds_names, paste("fold", i, sep = "")) 
  #   }
  # }
  return(folds)
}

cross_validation <- function(strat_data, k = 3) {
  results <- list()
  for (i in names(strat_data) ) {
    train_data <- strat_data[[i]]
    train_matrix <- as.matrix(train_data) #matrix of validation
    validation_data <- strat_data
    validation_data$i <- NULL #eliminates the validation dataset from the list of all folds
    validation_matrix <- as.matrix(unlist(validation_data, recursive = FALSE))
    print(as.factor(colnames(train_data)))
    classifier <- SWAP.KTSP.Train(train_matrix, as.factor(colnames(train_data)), krange = k) #trains classifier
    prediction <- SWAP.KTSP.Classify(validation_matrix, classifier) #uses classifier to ... classify the validation set
    results <- append(results, prediction)
  }
}

error_metric = function(table) {
  TP =table[1,1]
  TN =table[2,2]
  FP =table[1,2]
  FN =table[2,1]
  precision = (TP) / (TP + FP)
  accuracy_model  = (TP + TN) / (TP + TN + FP + FN)
  recall = (TP) / (TP + FN)
  print( paste("Precision value of the model: ", round(precision, 3)))
  print( paste("Accuracy of the model: ", round(accuracy_model, 3)))
  print( paste("Recall of the model: ", round(recall, 3)))
}


data_by_class = sep_by_class(df, axis = 2) #Marche pas avec axis = 1

data_by_class_by_fold <- strat_kfold(df, k = 3)

cross_validation(data_by_class_by_fold)

df_T <- t(df)

data_by_class <- sep_by_class(df, axis = 2)

data_by_class_T <- sep_by_class(df_T, axis = 1) #marche pas encore

## stratified k-fold must be the way to validate the model, because of the small dataset and imbalance between the 2 classes