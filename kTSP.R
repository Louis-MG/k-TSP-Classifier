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
col_names = gsub('^SLS.*', 'autre', col_names) #uses regex to change lines starting with SLS (followed by any character any number of times), by irritant
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


############################################
#
#   CROSS-VALIDATION 
#
############################################


### KFOLD STRATIFIE

#la premiere fonction separe le jeu de donnees en dataframes par categories, en calculant la proportion du jeu de donnees initiale pour chaque categorie

sep_by_label = function(df) {
  #
  #' df  : dataframe. Genes as rows and individuals as columns.
  #' This function separates the dataframe into a list of dataframes, each list corresponding to a class
  #
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  list_prop <- list() #list of the proportion of each class
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
  for (i in S3_object$labels) {
   colnames(S3_object[[i]]) <- rep(i, dim(S3_object[[i]])[2])
  }
  return(S3_object) #returns S3 class object with proportions for each class, dataframes for the data from each class, in the same order 
}

#WORKS GREAT HELL TEAH BABY
rename_ <- function(x, i) {
  #
  #' x : dataframe
  #' i : a class label
  #' This function renames the classes correctly after some manipulation leads to the addition of unwanted characters to column's names
  #
  col_names <- colnames(x) 
  col_names <- gsub(paste('^.*', i[1], '.*$', sep = ""), i[1], col_names)
  col_names <- gsub(paste('^.*', i[2], '.*$', sep = ""), i[2], col_names)
  colnames(x) <- col_names
  return(x)
}

strat_kfold = function(df, k) {
  #
  #' df : a dataframe
  #' k : numnber of subsets wanted
  #' This function uses the data separated by class to yield a list of dataframes, each being a sub-dataset with the original proportions of each class. 
  #
  if (is.data.frame(df) == FALSE) {
    stop("ERROR: df must be a dataframe.")
  }
  data_by_label <- sep_by_label(df) #splits dataframe by class in sub dataframe as S3 class object attributes
  sets_index = list() #list of indexes for the k subsets for kfold cv
  sets_index <- lapply(data_by_label[3:length(data_by_label)], function(x) split(sample(c(1:dim(x)[2]), dim(x)[2]), sort(c(1:dim(x)[2])%%k +1))) #list, for each class there are k lists containing the indexes 
  names(sets_index) <- as.list(data_by_label$labels)
  data_by_fold <- lapply(data_by_label$labels, function(x) lapply(sets_index[[x]], function(y) data_by_label[[x]][,y])) #here the indexes are used to make a list of lists : level one is the class, level 2 (sublevel) is the fold
  names(data_by_fold) <- as.list(data_by_label$labels)
  folds <- list() #initialises the final object
  folds_names <- list()
  for (i in c(1:k)) {
    folds <- append(folds, list(as.data.frame(lapply(data_by_label$labels, function(x) data_by_fold[[x]][i]) )) ) #adds class1 and class2 datasets together to for k sub datasets with equal proportions of each class
    folds <- lapply(folds, function(x) rename_(x, data_by_label$labels)) 
    folds_names <- append(folds_names, paste("fold", i, sep = "")) 
  }
  names(folds) <- folds_names
  return(folds)
}


classify <- function(df_by_fold, fold, k = 3, N) {
  #
  #' df_by_fold : list of dataframes. Each dataframe is a fold for k-fold cross-validation
  #' fold : element of the list df_by_fold. Stratified subset of data to use as validation data for the cross-validation round.
  #' k : number of folds for cross-validation
  #' N : number of TSPs to retain.
  #' This function uses the fold as validation data and the other subsets as training data. Returns a confusion matrix, the metrics and the k TSPs from the classifier.
  #
  validation_data <- fold #validation set
  validation_matrix <- as.matrix(validation_data)
  reality <- as.factor(colnames(validation_data)) #truth about validation set labels
  train_data <- df_by_fold
  train_data$i <- NULL #eliminates the validation dataset from the list of all folds
  train_data <- do.call("cbind", train_data) #binds together the training folds
  train_data <- rename_(train_data, data_by_label$labels)
  train_matrix <- as.matrix(train_data)
  classifier <- SWAP.KTSP.Train(train_matrix, as.factor(colnames(train_data)), krange = N) #trains classifier
  prediction <- SWAP.KTSP.Classify(validation_matrix, classifier) #uses classifier to ... classify the validation set
  result_table <- table(prediction, reality)
  print(result_table)
  error_metrics(result_table)
  print(classifier$TSPs)
}


cross_validation <- function(df, k = 3, N = 3) {
  #
  #' df : dataframe
  #' k : number of folds for cross-validation
  #' N : number of TSPs to retain 
  #' This function applies cross-validation using k-TSP to the data. It returns the confusion matrix, metrics and N TSPs for each cross-validation round.
  #
  data_by_label <- sep_by_label(df) 
  data_by_fold <- strat_kfold(df, k)
  lapply(data_by_fold, function(x) classify(data_by_fold, x, k, N))
}

error_metrics = function(table) {
  #
  #' table : confusion matrix (actually a table)
  #' This function calculates precision, accuracy, and recall scores of a classification with known labels/classes.
  #
  TP <- table[1,1]
  TN <- table[2,2]
  FP <- table[1,2]
  FN <- table[2,1]
  precision <- (TP) / (TP + FP)
  accuracy_model  <- (TP + TN) / (TP + TN + FP + FN)
  recall <- (TP) / (TP + FN)
  print( paste("Precision value of the model: ", round(precision, 3)))
  print( paste("Accuracy of the model: ", round(accuracy_model, 3)))
  print( paste("Recall of the model: ", round(recall, 3)))
}


data_by_label = sep_by_label(df) #Marche pas avec axis = 1

data_by_fold <- strat_kfold(df, k = 3)

cross_validation(df, k = 4, N = 5)

############################################
#
#   CLASSIFIER
#
############################################

compare_rules <- function(individual, confidence = FALSE) {
  #
  #' individual : columns name of a dataframe that we want to assign to a class
  #' confidence : boolean. Shows confidence in class attribution
  #' This function uses kTSP rules from the stratified cross-validation to assign individuals to a class of eczema.
  #
  autre = 0
  allergique = 0
  if (df["GBP1", individual] < df["LY6K", individual]) {
    autre = autre + 1
  } else {
    allergique = allergique + 1
  }
  if (df["PLSCR1.1", individual] < df["SCCPDH", individual]) {
    autre = autre + 1
  } else {
    allergique = allergique + 1
  }
  if (df["CDKN2D", individual] < df["PLA2G4F", individual]) {
    autre = autre + 1
  } else {
    allergique = allergique + 1
  }
  if (df["DTX3L.1", individual] < df["ANKK1", individual]) {
    autre = autre + 1
  } else {
    allergique = allergique + 1
  }
  if (df["ERAP2.1", individual] < df["GSTA3", individual]) {
    autre = autre + 1
  } else {
    allergique = allergique + 1
  }
  if (autre > allergique) {
    if (confidence) {
      print(autre/5)
      return("autre")
    } else {
      return("autre")
    }
  } else {
    if (confidence) {
      print(allergique/5)
      return("allergique")
    } else {
      return("allergique")
    }
  }
}

eczema_classifier <- function(df, confidence = FALSE) {
  #
  #' df : dataframe. Genes must be in rows and individuals in columns
  #' confidence : boolean. Shows confidence in class attribution
  #' This function is the actual classifier to distinguish allergic eczema from other types. Returns a list with classes for each individual.
  #
  df <- df[ which(rownames(df) %in% c("LY6K", "GBP1", "SCCPDH", "PLSCR1.1", "PLA2G4F", "CDKN2D", "ANKK1", "DTX3L.1", "GSTA3", "ERAP2.1")) , ]
  classes <- lapply(colnames(df), function(x) compare_rules(x, confidence))
  return(classes)
}

test = eczema_classifier(df)
test = eczema_classifier(df, confidence = TRUE) #permet de visualiser la confiance en l'attribution de la classe. Compris dans [0;1]
test == colnames(df) #permet de comparer les prédictions avec les classes déjà attribuées du jeu de données
sum(test == colnames(df)) #la valeur TRUE vaut 1, et FALSE vaut 0, donc faire la sum d'un array de TRUE et FALSE permet de connaître le nombre de TRUE
sum(test == colnames(df))/dim(df)[2]#par extension, on peut transformer cette somme en % pour connaître la justesse d'autant de prédicitons que l'on veut (sans regarder "à la main", ce qui est horrible sur bcp de données)
