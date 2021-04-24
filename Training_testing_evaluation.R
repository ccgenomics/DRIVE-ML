########################################################################################################################
### LOAD PACKAGES REQUIRED FOR DATA PRE-PROCESSING AND MACHINE LEARNING 
########################################################################################################################

### DATA MANIPULATION PACKAGES
library(plyr) 
library(tidyverse)

### MACHINE LEARNING/DEEP LEARNING PACKAGES
library(caret) #required for training ML models
library(randomForest) #required for random forest
library(multcomp) #required for caret package
library(party) #required for caret package
library(doParallel) #required for parallel processing
library(e1071) #required for caret package
library(pROC) #required for computing the AUC-ROC
library(rpart.plot) #required for decision trees
library(ROCR) #required for plotting AUC-ROC
library(xgboost) #required for extreme gradient boosting
library(gbm) #required for stochastic gradient boosting
library(caTools) #required for logistic regression
library(MLeval) #required for ML model evaluation
library(mltest) #required for ML evaluation based on confusion matrix
library(kknn) #required for knn
library(kernlab) #required for SVM
library(DMwR) #required for SMOTE (oversampling/undersampling if needed)
library(UBL) #required for undersampling
library(ranger) #required for random forest
library(lattice)
library(tuneRanger)
library(ggROC)
library(base)
library(glmnet) #required for logistic regression
library(h2o) #required for connection to a Java VM for parallel processing (version older than 13 required - 11 used here)
library(xgboost) #required for extreme gradient boosting
library(varImp) #required for variable importance
library(hmeasure) #required for confusion matrix of CHASMplus probabilities for performance metrics comparison during benchmarking
library(RSNNS) #required for MLP
library(keras) #required for MLP - other versions

### VISUALISATION PACKAGES
library(ggplot2) #required for graphics
library(ggpubr) #required for merging feature importance figures





########################################################################################################################
### LOAD TRAINING DATA (GENIE) AND PRE-PROCESS FOR A ML FORMAT
########################################################################################################################

#load full dataset
data_GENIE <- read.table("data_ready_for_ML_GENIE.txt", header=TRUE, sep="")

#remove rows with missing data
data_GENIE <- data_GENIE[complete.cases(d), ]

#copy the row names as first column
library(data.table)
setDT(data_GENIE, keep.rownames = TRUE)[]

#if needed, the data can be reduced for testing on a small subset of data before pushing a job to the VM
#library(dplyr)
#reduced_GENIE <- d %>%
#group_by(labels_driver_mutation) %>%
#sample_frac(186060/nrow(data_GENIE))

#perform undersampling for class imbalance so that no of drivers is equal to the no of passengers
under_samp_reduced_GENIE <- RandUnderClassif(labels_driver_mutation~., as.data.frame(data_GENIE), C.perc="balance")
summary(under_samp_reduced_GENIE)

#make the first column as row names and then delete the first column
row.names(under_samp_reduced_GENIE) <- under_samp_reduced_GENIE$rn
under_samp_reduced_GENIE$rn <- NULL

#create a new file for factor level processing
factor_dataset_GENIE <- under_samp_reduced_GENIE

#create factor columns for each categorical variable and labels
cols_to_be_factorised_GENIE <- c("labels_driver_mutation",
                                 "cellular.component",
                                 "development",
                                 "DNA.damage",
                                 "immune",
                                 "metabolic",
                                 "pathway",
                                 "proliferation",
                                 "signaling",
                                 "NA_hallmarks",
                                 "acetylation",
                                 "caspase.cleavage",
                                 "di.methylation",
                                 "methylation",
                                 "mono.methylation",
                                 "N.Glycosylation",
                                 "O.GalNAc",
                                 "O.GlcNAc",
                                 "phosphorylation",
                                 "succinylation",
                                 "sumoylation",
                                 "tri.methylation",
                                 "ubiquitylation")
factor_dataset_GENIE[cols_to_be_factorised_GENIE] <- lapply(factor_dataset_GENIE[cols_to_be_factorised_GENIE], as.factor)

#rename the dataset back with the same name
under_samp_reduced_GENIE <- factor_dataset_GENIE





########################################################################################################################
### LOAD TESTING DATA (TCGA) AND PRE-PROCESS FOR A ML FORMAT
########################################################################################################################

#load data
data_TCGA <- read.csv("data_ready_for_ML_TCGA.txt", header=TRUE, sep="")

#remove rows with missing data
data_TCGA <- data_TCGA[complete.cases(data_TCGA), ]

#copy the row names as first column
library(data.table)
setDT(data_TCGA, keep.rownames = TRUE)[]

#if needed, the data can be reduced for testing on a small subset of data before pushing a job to the VM
#library(dplyr)
#reduced_TCGA <- d %>%
#group_by(labels_driver_mutation) %>%
#sample_frac(186060/nrow(data_TCGA))

#perform undersampling for class imbalance so that no of drivers is equal to the no of passengers
under_samp_reduced_TCGA <- RandUnderClassif(labels_driver_mutation~., as.data.frame(data_TCGA), C.perc="balance")
summary(under_samp_reduced_TCGA)

#make the first column as row names and then delete the first column
row.names(under_samp_reduced_TCGA) <- under_samp_reduced_TCGA$rn
under_samp_reduced_TCGA$rn <- NULL


#create a new file for factor level processing
factor_dataset_TCGA <- under_samp_reduced_TCGA

#create factor columns for each categorical variable and labels
cols_to_be_factorised_TCGA <- c("labels_driver_mutation",
                                "cellular.component",
                                "development",
                                "DNA.damage",
                                "immune",
                                "metabolic",
                                "pathway",
                                "proliferation",
                                "signaling",
                                "NA_hallmarks",
                                "acetylation",
                                "caspase.cleavage",
                                "di.methylation",
                                "methylation",
                                "mono.methylation",
                                "N.Glycosylation",
                                "O.GalNAc",
                                "O.GlcNAc",
                                "phosphorylation",
                                "succinylation",
                                "sumoylation",
                                "tri.methylation",
                                "ubiquitylation")
factor_dataset_TCGA[cols_to_be_factorised_TCGA] <- lapply(factor_dataset_TCGA[cols_to_be_factorised_TCGA], as.factor)

#rename the dataset back with the same name
under_samp_reduced_TCGA <- factor_dataset_TCGA




########################################################################################################################
### LOAD BENCHMARK DATASET (MutaGene) AND PRE-PROCESS FOR A ML FORMAT
########################################################################################################################

#load full dataset
data_BENCHMARK <- read.table("data_ready_for_ML_BENCHMARK.txt", header=TRUE, sep="\t")

#remove rows with missing data
data_BENCHMARK <- data_BENCHMARK[complete.cases(data_BENCHMARK), ]

#copy the row names as first column
library(data.table)
setDT(data_BENCHMARK, keep.rownames = TRUE)[]

#if needed, the data can be reduced for testing on a small subset of data before pushing a job to the VM
#library(dplyr)
#reduced_BENCHMARK <- d %>%
#group_by(labels_driver_mutation) %>%
#sample_frac(186060/nrow(reduced_BENCHMARK))

#perform undersampling for class imbalance so that no of drivers is equal to the no of passengers
under_samp_reduced_BENCHMARK <- RandUnderClassif(Class_label_MutaGene~., as.data.frame(data_BENCHMARK), C.perc="balance")
summary(under_samp_reduced_BENCHMARK)

#make the first column as row names and then delete the first column
row.names(under_samp_reduced_BENCHMARK) <- make.names(under_samp_reduced_BENCHMARK$VEP_merge_col2, unique=TRUE)
under_samp_reduced_BENCHMARK$rn <- NULL
under_samp_reduced_BENCHMARK$VEP_merge_col2 <- NULL

#create a new file for factor level processing
factor_dataset_BENCHMARK <- under_samp_reduced_BENCHMARK

#create factor columns for each categorical variable and labels
cols_to_be_factorised_BENCHMARK <- c("Class_label_MutaGene",
                                     "cellular.component",
                                     "development",
                                     "DNA.damage",
                                     "immune",
                                     "metabolic",
                                     "pathway",
                                     "proliferation",
                                     "signaling",
                                     "NA_hallmarks",
                                     "acetylation",
                                     "caspase.cleavage",
                                     "di.methylation",
                                     "methylation",
                                     "mono.methylation",
                                     "N.Glycosylation",
                                     "O.GalNAc",
                                     "O.GlcNAc",
                                     "phosphorylation",
                                     "succinylation",
                                     "sumoylation",
                                     "tri.methylation",
                                     "ubiquitylation")
factor_dataset_BENCHMARK[cols_to_be_factorised_BENCHMARK] <- lapply(factor_dataset_BENCHMARK[cols_to_be_factorised_BENCHMARK], as.factor)

#rename the dataset back with the same name
under_samp_reduced_BENCHMARK <- factor_dataset_BENCHMARK





########################################################################################################################
### PREPARE TRAINING (GENIE), TEST (TCGA) AND BENCHMARK (MutaGene) DATA FOR A DL FORMAT - necessary for mlp models
########################################################################################################################

#dummify GENIE data (for one-hot encoding)
dummy_data_GENIE <- dummyVars("~ .", data=under_samp_reduced_GENIE[,1:51])
transformed_dum_data_GENIE <- data.frame(predict(dummy_data_GENIE, newdata=under_samp_reduced_GENIE[,1:51]))
transformed_dum_data_GENIE$labels_driver_mutation <- under_samp_reduced_GENIE$labels_driver_mutation

#dummify TCGA data (for one-hot encoding)
dummy_data_TCGA <- dummyVars("~ .", data=under_samp_reduced_TCGA[,1:51])
transformed_dum_data_TCGA <- data.frame(predict(dummy_data_TCGA, newdata=under_samp_reduced_TCGA[,1:51]))
transformed_dum_data_TCGA$labels_driver_mutation <- under_samp_reduced_TCGA$labels_driver_mutation


#dummify BENCHMARK data (for one-hot encoding) - ONLY NEEDED IF BENCHMARKING USING mlp model
dummy_data_BENCHMARK <- dummyVars("~ .", data=under_samp_reduced_BENCHMARK[,1:51])
transformed_dum_data_BENCHMARK <- data.frame(predict(dummy_data_BENCHMARK, newdata=under_samp_reduced_BENCHMARK[,1:51]))
transformed_dum_data_BENCHMARK$labels_driver_mutation <- under_samp_reduced_BENCHMARK$labels_driver_mutation





########################################################################################################################
### RANDOM FOREST
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation (method of resampling)
                              number=10, #k=10
                              savePredictions="final", #save predictions for the optimal tuning parameters
                              summaryFunction=twoClassSummary, #compute summary metrics across the resamples
                              classProbs=TRUE,
                              verboseIter=TRUE) #add class probabilities & predicted values for the classification model in each resample

model_rf <- caret::train(x=under_samp_reduced_GENIE[,1:51], #predictors
                         y=under_samp_reduced_GENIE[,52], #outcome
                         method="ranger", #random forest
                         metric="ROC", #metric used for selecting the optimal model
                         tuneLength=10, #how many levels of parameter tuning should be tried
                         preProcess=c("scale"), #pre-processing of predictors using scale method (divide the values by their st dev)
                         trControl=train_control,
                         importance="impurity") #define how the train function should work 

#extract the final model parameters
model_rf$bestTune

#performance metrics for prediction average from all folds 
cm_rf <- caret::confusionMatrix(model_rf$pred[order(model_rf$pred$rowIndex),4], #order the 4th column of the list (i.e. predictions on driver mutations)
                                under_samp_reduced_GENIE$labels_driver_mutation, #original GENIE data
                                mode="everything") #extract all performance measurements
cm_rf



### TESTING &  PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
rf_classes <- predict(model_rf, newdata=under_samp_reduced_TCGA[,1:51])

#predict using raw labels and prediction probabilities
rf_predicted_labels <- as.data.frame(predict(model_rf, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
rf_predicted_prob <- predict(model_rf, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

#extract performance metrics of the prediction
rf_cm <- caret::confusionMatrix(data=rf_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
rf_cm

#extract AUC of the prediction
rf_ROC_pred <- ROCR::prediction(rf_predicted_prob$driver,
                                under_samp_reduced_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
rf_performance <- ROCR::performance(rf_ROC_pred, "auc")
rf_performance@y.values[[1]]



### BENCHMARKING
#using the model, predict on the BENCHMARK dataset for driver vs. passenger
rf_classes_benchmark <- predict(model_rf, newdata=under_samp_reduced_BENCHMARK[,1:51])

#predict using raw labels and prediction probabilities
rf_predicted_labels_benchmark <- as.data.frame(predict(model_rf, newdata=under_samp_reduced_BENCHMARK[,1:51], type="raw"))
rf_predicted_prob_benchmark <- predict(model_rf, newdata=under_samp_reduced_BENCHMARK[,1:51], type="prob")

#extract performance metrics of the prediction
rf_cm_benchmark <- caret::confusionMatrix(data=rf_classes_benchmark, under_samp_reduced_BENCHMARK$Class_label_MutaGene, mode="everything")
rf_cm_benchmark

#extract AUC of the prediction
rf_ROC_pred_benchmark <- ROCR::prediction(rf_predicted_prob_benchmark$driver,
                                          under_samp_reduced_BENCHMARK$Class_label_MutaGene, label.ordering = c("passenger", "driver"))
rf_performance_benchmark <- ROCR::performance(rf_ROC_pred_benchmark, "auc")
rf_performance_benchmark@y.values[[1]]

#extract performance metrics of CHASMplus probabilities
CHASMplus_performance_metrics <- HMeasure(under_samp_reduced_BENCHMARK$Class_label_MutaGene, under_samp_reduced_BENCHMARK$CHASMplus_score)
summary(CHASMplus_performance_metrics, show.all=TRUE)



### train another time with importance set to permutation/mean decrease accuracy
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation (method of resampling)
                              number=10, #k=10
                              savePredictions="final", #save predictions for the optimal tuning parameters
                              summaryFunction=twoClassSummary, #compute summary metrics across the resamples
                              classProbs=TRUE,
                              verboseIter=TRUE) #add class probabilities & predicted values for the classification model in each resample

model_rf2 <- caret::train(x=under_samp_reduced_GENIE[,1:51], #predictors
                          y=under_samp_reduced_GENIE[,52], #outcome
                          method="ranger", #random forest
                          metric="ROC", #metric used for selecting the optimal model
                          tuneLength=10, #how many levels of parameter tuning should be tried
                          preProcess=c("scale"), #pre-processing of predictors using scale method (divide the values by their st dev)
                          trControl=train_control,
                          importance="permutation") #define how the train function should work 


### FEATURE IMPORTANCE based on impurity/mean decrease gini
rf_var_imp_gini <- caret::varImp(model_rf) #calculate feature importance based in impurity
rf_var_imp_gini_df <- as.data.frame(rf_var_imp_gini$importance) #store importance results in a data frame
rf_var_imp_gini_df$feature_name <- rownames(rf_var_imp_gini_df) #make rownames as column
rownames(rf_var_imp_gini_df) <- NULL #delete rownames

#create a vector with gene level features for category
vector_gene_level_features <- c("max_no_interactions", "cellular.component", 
                                "development",                          "DNA.damage",                          
                                "immune",                               "metabolic",                           
                                "pathway",                              "proliferation" ,                      
                                "signaling",                            "NA_hallmarks",                        
                                "acetylation",                          "caspase.cleavage" ,                   
                                "di.methylation",                       "methylation" ,                        
                                "mono.methylation" ,                    "N.Glycosylation",                     
                                "O.GalNAc" ,                            "O.GlcNAc" ,                           
                                "phosphorylation",                      "succinylation" ,                      
                                "sumoylation",                          "tri.methylation" ,                    
                                "ubiquitylation",                       "gene.length" ,                        
                                "silent" ,                              "nonsense" ,                           
                                "splice.site",                          "missense" ,                           
                                "recurrent.missense" ,                  "normalized.missense.position.entropy",
                                "frameshift.indel",                     "inframe.indel",                       
                                "normalized.mutation.entropy",          "Mean.Missense.MGAEntropy" ,           
                                "Mean.VEST.Score" ,                     "lost.start.and.stop",                 
                                "missense.to.silent",                   "non.silent.to.silent" ,               
                                "expression_CCLE" ,                     "replication_time" ,                   
                                "HiC_compartment",                      "gene_betweeness"  ,                   
                                "gene_degree"  ,                        "oncogene.score",                      
                                "tsg.score" ,                           "other.score",                         
                                "driver.score")

#create a vector with mutation level features for category
vector_mutation_level_features <- c("SIFT",                                
                                    "PolyPhen", "Condel",                              
                                    "average_rank")
#add the feature category
rf_var_imp_gini_df$feature_category <- ifelse(rf_var_imp_gini_df$feature_name %in% vector_gene_level_features, "Gene level", "Mutation level")

#sort the dataframe by the impurity overall
rf_var_imp_gini_df_ordered <- rf_var_imp_gini_df[order(-rf_var_imp_gini_df$Overall),]

#plot the variable importance for impurity
library(ggplot2)
bold.text <- element_text(face="bold") #define a bold text for adding to the x and y axes


feature_imp_MDI <- ggplot(rf_var_imp_gini_df_ordered, 
                          frame.plot=FALSE, #remove plot frame
                          scale_y_discrete(expand=c(0,0)), #set values for discrete y-axis
                          mapping = aes(x=reorder(feature_name, Overall), #set figure aesthetics
                                        y=Overall, 
                                        color=(as.factor(feature_category)))) +
  geom_point()+
  geom_segment(aes(x=feature_name, xend=feature_name, y=0, yend=Overall)) +
  scale_color_discrete(name="Feature level") + #color by feature level
  ylab("Importance based on impurity") + #add Importance based on impurity label
  xlab("Feature name") + #add Feature name labels
  theme(axis.title = bold.text) + #bold axes labels
  coord_flip() #flip coordinates so that x-axis becomes y-axis and vice-versa




### FEATURE IMPORTANCE based on mean decrease accuracy
rf_var_imp_acc <- caret::varImp(model_rf2) #calculate feature importance based in impurity
rf_var_imp_acc_df <- as.data.frame(rf_var_imp_acc$importance) #store importance results in a data frame
rf_var_imp_acc_df$feature_name <- rownames(rf_var_imp_acc_df) #make rownames as column
rownames(rf_var_imp_acc_df) <- NULL #delete rownames

#create a vector with gene level features for category
vector_gene_level_features <- c("max_no_interactions", "cellular.component", 
                                "development",                          "DNA.damage",                          
                                "immune",                               "metabolic",                           
                                "pathway",                              "proliferation" ,                      
                                "signaling",                            "NA_hallmarks",                        
                                "acetylation",                          "caspase.cleavage" ,                   
                                "di.methylation",                       "methylation" ,                        
                                "mono.methylation" ,                    "N.Glycosylation",                     
                                "O.GalNAc" ,                            "O.GlcNAc" ,                           
                                "phosphorylation",                      "succinylation" ,                      
                                "sumoylation",                          "tri.methylation" ,                    
                                "ubiquitylation",                       "gene.length" ,                        
                                "silent" ,                              "nonsense" ,                           
                                "splice.site",                          "missense" ,                           
                                "recurrent.missense" ,                  "normalized.missense.position.entropy",
                                "frameshift.indel",                     "inframe.indel",                       
                                "normalized.mutation.entropy",          "Mean.Missense.MGAEntropy" ,           
                                "Mean.VEST.Score" ,                     "lost.start.and.stop",                 
                                "missense.to.silent",                   "non.silent.to.silent" ,               
                                "expression_CCLE" ,                     "replication_time" ,                   
                                "HiC_compartment",                      "gene_betweeness"  ,                   
                                "gene_degree"  ,                        "oncogene.score",                      
                                "tsg.score" ,                           "other.score",                         
                                "driver.score")

#create a vector with mutation level features for category
vector_mutation_level_features <- c("SIFT",                                
                                    "PolyPhen", "Condel",                              
                                    "average_rank")
#add the feature category
rf_var_imp_acc_df$feature_category <- ifelse(rf_var_imp_acc_df$feature_name %in% vector_gene_level_features, "Gene level", "Mutation level")

#sort the dataframe by the permutations overall
rf_var_imp_acc_df_ordered <- rf_var_imp_acc_df[order(-rf_var_imp_acc_df$Overall),]

#plot the variable importance for permutations
library(ggplot2)
bold.text <- element_text(face="bold") #define a bold text for adding to the x and y axes
feature_imp_MDA <- ggplot(rf_var_imp_acc_df_ordered, 
                          frame.plot=FALSE, #remove plot frame
                          mapping = aes(x=reorder(feature_name, Overall), #set figure aesthetics
                                        y = Overall, 
                                        color = (as.factor(feature_category)))) + 
  geom_point()+
  geom_segment(aes(x=feature_name, xend=feature_name, y=0, yend=Overall)) +
  scale_color_discrete(name="Feature level") + #color by feature level
  ylab("Importance based on permutations") + #add Importance based on permutations/accuracy label
  xlab("") + #remove Feature name from the axis name
  theme(axis.title = bold.text) + #bold axes labels
  coord_flip() #flip coordinates so that x-axis becomes y-axis and vice-versa

#MERGE BOTH FEATURE IMPORTANCE FIGURE PARTS (MDI & MDA)
tiff("RF_importance.tiff", width = 3200, height = 1900, units = 'px', res = 300)

ggpubr::ggarrange(feature_imp_MDI, feature_imp_MDA, #merge both parts of random forest feature importance
                  labels=c("a)", "b)"), #add labels for each part
                  ncol=2, nrow=1, #specity no of rows and columns
                  common.legend = TRUE, #make just one common legend
                  legend="right") #set position of the legend as on the right

dev.off()





########################################################################################################################
#### DECISION TREE
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation
                              number=10,
                              savePredictions="final",
                              summaryFunction=twoClassSummary,
                              classProbs=TRUE,
                              verboseIter=TRUE)

model_dt <- caret::train(x=under_samp_reduced_GENIE[,1:51],
                         y=under_samp_reduced_GENIE[,52],
                         method="rpart",
                         metric="ROC",
                         tuneLength=10,
                         preProcess=c("scale"),
                         trControl=train_control)

#extract the final model parameters
model_dt$bestTune


#performance metrics for prediction average from all folds 
cm_dt <- caret::confusionMatrix(model_dt$pred[order(model_dt$pred$rowIndex),2], #order the 2nd column of the list (i.e. predictions on driver mutations)
                                under_samp_reduced_GENIE$labels_driver_mutation, #original GENIE data
                                mode="everything")
cm_dt



#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
dt_classes <- predict(model_dt, newdata=under_samp_reduced_TCGA[,1:51])

#predict using raw labels and prediction probabilities
dt_predicted_labels <- as.data.frame(predict(model_dt, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
dt_predicted_prob <- predict(model_dt, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

#extract performance metrics of the prediction
dt_cm <- caret::confusionMatrix(data=dt_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
dt_cm

#extract AUC of the prediction
dt_ROC_pred <- ROCR::prediction(dt_predicted_prob$driver,
                                under_samp_reduced_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
dt_performance <- ROCR::performance(dt_ROC_pred, "auc")
dt_performance@y.values[[1]]





########################################################################################################################
#### EXTREME GRADIENT BOOSTING
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation
                              number=10,
                              savePredictions=TRUE,
                              summaryFunction=twoClassSummary,
                              classProbs=TRUE,
                              verboseIter=TRUE)


model_egb <- caret::train(labels_driver_mutation ~ .,
                          data=under_samp_reduced_GENIE,
                          method="xgbTree",
                          metric="ROC",
                          tuneLength=10,
                          nthread=1,
                          preProcess=c("scale"),
                          trControl=trainControl(method="cv", #k-fold cross validation
                                                 number=10,
                                                 savePredictions=TRUE,
                                                 summaryFunction=twoClassSummary,
                                                 classProbs=TRUE, 
                                                 verboseIter=TRUE))

#load the model (if it was produced and downloaded from a VM)
#load("model_egb.RData")


#extract the final model parameters
model_egb$bestTune


#performance metrics for prediction average from all folds
cm_egb <- caret::confusionMatrix(model_egb$pred[order(model_egb$pred$rowIndex),1][1:24594], #order the 1ST column of the list (i.e. predictions on driver mutations)
                                 under_samp_reduced_GENIE$labels_driver_mutation, 
                                 mode="everything")
cm_egb

#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
egb_classes <- predict(model_egb, newdata=under_samp_reduced_TCGA[,1:51])

#predict using raw labels and prediction probabilities
egb_predicted_labels <- as.data.frame(predict(model_egb, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
egb_predicted_prob <- predict(model_egb, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

#extract performance metrics of the prediction
egb_cm <- caret::confusionMatrix(data=egb_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
egb_cm

#extract AUC of the prediction
egb_ROC_pred <- ROCR::prediction(egb_predicted_prob$driver,
                                 under_samp_reduced_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
egb_performance <- ROCR::performance(egb_ROC_pred, "auc")
egb_performance@y.values[[1]]




########################################################################################################################
### LOGISTIC REGRESSION
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
model_lr <- caret::train(labels_driver_mutation ~ .,
                         data=under_samp_reduced_GENIE,
                         method="glmnet",
                         metric="ROC",
                         tuneLength=10,
                         preProcess=c("scale"),
                         trControl=trainControl(method="cv", #k-fold cross validation
                                                number=10,
                                                savePredictions=TRUE,
                                                summaryFunction=twoClassSummary,
                                                classProbs=TRUE,
                                                verboseIter=TRUE))

#load the model (if it was previously run on a VM)
#load("model_lr.RData")

#extract the final model parameters
model_lr$bestTune

cm_lr <- caret::confusionMatrix(model_lr$pred[order(model_lr$pred$rowIndex),1][1:24594], #order the 1st column of the list (i.e. predictions on driver mutations)
                                under_samp_reduced_GENIE$labels_driver_mutation, 
                                mode="everything")
cm_lr

#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
lr_classes <- predict(model_lr, newdata=under_samp_reduced_TCGA[,1:51])

#predict using raw labels and prediction probabilities
lr_predicted_labels <- as.data.frame(predict(model_lr, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
lr_predicted_prob <- predict(model_lr, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

#extract performance metrics of the prediction
lr_cm <- caret::confusionMatrix(data=lr_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
lr_cm

#extract AUC of the prediction
lr_ROC_pred <- ROCR::prediction(lr_predicted_prob$driver,
                                under_samp_reduced_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
lr_performance <- ROCR::performance(lr_ROC_pred, "auc")
lr_performance@y.values[[1]]





########################################################################################################################
### SVM
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation
                              number=10,
                              savePredictions=TRUE,
                              summaryFunction=twoClassSummary,
                              classProbs=TRUE,
                              verboseIter=TRUE)


model_svm <- caret::train(labels_driver_mutation ~ .,
                          data=under_samp_reduced_GENIE,
                          method="svmRadial",
                          metric="ROC",
                          tuneLength=5,
                          preProcess=c("scale"),
                          trControl=train_control,
                          threshold=0.3)

#extract the final model parameters
model_svm$bestTune

cm_svm <- caret::confusionMatrix(model_svm$pred[order(model_svm$pred$rowIndex),1][1:24594], 
                                 under_samp_reduced_GENIE$labels_driver_mutation, 
                                 mode="everything")
cm_svm

#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
svm_classes <- predict(model_svm, newdata=under_samp_reduced_TCGA[,1:51])

#predict using raw labels and prediction probabilities
svm_predicted_labels <- as.data.frame(predict(model_svm, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
svm_predicted_prob <- predict(model_svm, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

#extract performance metrics of the prediction
svm_cm <- caret::confusionMatrix(data=svm_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
svm_cm

#extract AUC of the prediction
svm_ROC_pred <- ROCR::prediction(svm_predicted_prob$driver,
                                 under_samp_reduced_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
svm_performance <- ROCR::performance(svm_ROC_pred, "auc")
svm_performance@y.values[[1]]





########################################################################################################################
### KNN
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation
                              number=10,
                              savePredictions=TRUE,
                              summaryFunction=twoClassSummary,
                              classProbs=TRUE,
                              verboseIter=TRUE)


model_knn <- caret::train(labels_driver_mutation ~ .,
                          data=under_samp_reduced_GENIE,
                          method="knn",
                          metric="ROC",
                          tuneLength=10,
                          preProcess=c("scale"),
                          trControl=train_control)

#load the model (since it was generated using the VM)
#load("model_knn.RData")

#extract the final model parameters
model_knn$bestTune

#extract performance metrics of the training
cm_knn <- caret::confusionMatrix(model_knn$pred[order(model_knn$pred$rowIndex),1][1:24594], 
                                 under_samp_reduced_GENIE$labels_driver_mutation, 
                                 mode="everything")
cm_knn

#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
knn_classes <- predict(model_knn, newdata=under_samp_reduced_TCGA[,1:51])
levels(knn_classes)

#extract performance metrics of the prediction
knn_cm <- caret::confusionMatrix(data=knn_classes, under_samp_reduced_TCGA$labels_driver_mutation, mode="everything")
knn_cm


#predict using raw labels and prediction probabilities
knn_predicted_labels <- as.data.frame(predict(model_knn, newdata=under_samp_reduced_TCGA[,1:51], type="raw"))
knn_predicted_prob <- predict(model_knn, newdata=under_samp_reduced_TCGA[,1:51], type="prob")

under_samp_reduced_TCGA$numeric_labels <- as.numeric(under_samp_reduced_TCGA$labels_driver_mutation =="driver")

#extract AUC of the prediction
pred1 <- ROCR::prediction(knn_predicted_prob$driver, 
                          under_samp_reduced_TCGA$numeric_labels)
perf1 <- ROCR::performance(pred1, "tpr", "fpr")
plot(perf1)
auc <- ROCR::performance(pred1, "auc")
auc@y.values[[1]]





########################################################################################################################
### MULTILAYER PERCEPTRON (mlp simple)
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation (method of resampling)
                              number=10, #k=10
                              savePredictions="final", #save predictions for the optimal tuning parameters
                              summaryFunction=twoClassSummary, #compute summary metrics across the resamples
                              classProbs=TRUE,
                              verboseIter=TRUE) #add class probabilities & predicted values for the classification model in each resample

model_mlp <- caret::train(x=transformed_dum_data_GENIE[,1:73], #predictors
                          y=transformed_dum_data_GENIE[,74], #outcome
                          method="mlp", #multi-layer perceptron
                          metric="ROC", #metric used for selecting the optimal model
                          tuneLength=10, #how many levels of parameter tuning should be tried
                          preProcess=c("scale"), #pre-processing of predictors using scale method (divide the values by their st dev)
                          trControl=train_control) #define how the train function should work 

#extract the final model parameters
model_mlp$bestTune


#performance metrics for prediction average from all folds
cm_mlp <- caret::confusionMatrix(model_mlp$pred[order(model_mlp$pred$rowIndex),2], 
                                 transformed_dum_data_GENIE$labels_driver_mutation, 
                                 mode="everything")
cm_mlp

#TESTING & PERFORMANCE METRICS
#using the model, predict on the TCGA dataset for driver vs. passenger
mlp_classes <- predict(model_mlp, newdata=transformed_dum_data_TCGA[,1:73])

#predict using raw labels and prediction probabilities
mlp_predicted_labels <- as.data.frame(predict(model_mlp, newdata=transformed_dum_data_TCGA[,1:73], type="raw"))
mlp_predicted_prob <- predict(model_mlp, newdata=transformed_dum_data_TCGA[,1:73], type="prob")

#extract performance metrics of the prediction
mlp_cm <- caret::confusionMatrix(data=mlp_classes, transformed_dum_data_TCGA$labels_driver_mutation, mode="everything")
mlp_cm

#extract AUC of the prediction
mlp_ROC_pred <- ROCR::prediction(mlp_predicted_prob$driver,
                                 transformed_dum_data_TCGA$labels_driver_mutation, label.ordering = c("passenger", "driver"))
mlp_performance <- ROCR::performance(mlp_ROC_pred, "auc")
mlp_performance@y.values[[1]]





########################################################################################################################
### MULTILAYER PERCEPTRON (mlpKerasDropout)
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation (method of resampling)
                              number=10, #k=10
                              savePredictions="final", #save predictions for the optimal tuning parameters
                              summaryFunction=twoClassSummary, #compute summary metrics across the resamples
                              classProbs=TRUE,
                              verboseIter=TRUE) #add class probabilities & predicted values for the classification model in each resample

model_mlpKerasDropout <- caret::train(x=transformed_dum_data_GENIE[,1:73], #predictors
                                      y=transformed_dum_data_GENIE[,74], #outcome
                                      method="mlpKerasDropout", #multi-layer perceptron
                                      metric="ROC", #metric used for selecting the optimal model
                                      tuneLength=10, #how many levels of parameter tuning should be tried
                                      preProcess=c("scale"), #pre-processing of predictors using scale method (divide the values by their st dev)
                                      trControl=train_control) #define how the train function should work 


#extract the final model parameters
model_mlpKerasDropout$bestTune


#performance metrics for prediction average from all folds 
cm_mlpKerasDropout <- caret::confusionMatrix(model_mlpKerasDropout$pred[order(model_mlpKerasDropout$pred$rowIndex),2], 
                                             transformed_dum_data_GENIE$labels_driver_mutation, 
                                             mode="everything")
cm_mlpKerasDropout

#extract AUC-ROC
roc_mlp_cv <- pROC::roc(predictor = model_mlpKerasDropout$pred$driver, 
                        response = model_mlpKerasDropout$pred$obs,
                        levels =c("passenger", "driver"), #set levels for 1st) control vs. 2nd) case
                        direction="<",
                        plot=TRUE, 
                        col="orange", 
                        lwd=1, #line width of the ROC curve
                        legacy.axes=TRUE, #by default, roc function gives sensitivity vs. specificity, but we need 1-specificity instead
                        font=1, #make tick marks bolder 
                        las=1, #rotate y-axis labels to be 1=horizontal
                        grid=TRUE, #add a grid to the plot
                        xlab="False positive rate",
                        ylab="True positive rate",
                        font.lab=4 #font (bold the x-axis and y-axis)
                        )





########################################################################################################################
### MULTILAYER PERCEPTRON (mlpKerasDecayCost)
########################################################################################################################

### TRAINING & PERFORMANCE METRICS
set.seed(467)
train_control <- trainControl(method="cv", #k-fold cross validation (method of resampling)
                              number=10, #k=10
                              savePredictions="final", #save predictions for the optimal tuning parameters
                              summaryFunction=twoClassSummary, #compute summary metrics across the resamples
                              classProbs=TRUE,
                              verboseIter=TRUE) #add class probabilities & predicted values for the classification model in each resample

model_mlpKerasDecayCost <- caret::train(x=transformed_dum_data_GENIE[,1:73], #predictors
                                        y=transformed_dum_data_GENIE[,74], #outcome
                                        method="mlpKerasDecayCost", #multi-layer perceptron
                                        metric="ROC", #metric used for selecting the optimal model
                                        tuneLength=10, #how many levels of parameter tuning should be tried
                                        preProcess=c("scale"), #pre-processing of predictors using scale method (divide the values by their st dev)
                                        trControl=train_control) #define how the train function should work 

#extract the final model parameters
model_mlpKerasDecayCost$bestTune


#performance metrics for prediction average from all folds 
cm_mlpKerasDecayCost <- caret::confusionMatrix(model_mlpKerasDecayCost$pred[order(model_mlpKerasDecayCost$pred$rowIndex),2], 
                                               transformed_dum_data_GENIE$labels_driver_mutation, 
                                               mode="everything")
cm_mlpKerasDecayCost

#extract AUC-ROC
roc_mlp_cv <- pROC::roc(predictor = model_mlpKerasDecayCost$pred$driver, 
                        response = model_mlpKerasDecayCost$pred$obs,
                        levels =c("passenger", "driver"), #set levels for 1st) control vs. 2nd) case
                        direction="<",
                        plot=TRUE, 
                        col="orange", 
                        lwd=1, #line width of the ROC curve
                        legacy.axes=TRUE, #by default, roc function gives sensitivity vs. specificity, but we need 1-specificity instead
                        font=1, #make tick marks bolder 
                        las=1, #rotate y-axis labels to be 1=horizontal
                        grid=TRUE, #add a grid to the plot
                        xlab="False positive rate",
                        ylab="True positive rate",
                        font.lab=4 #font (bold the x-axis and y-axis)
                        )





########################################################################################################################
### EXTRACT AND PLOT AUC-ROC for TRAINING ON GENIE data (for all partitions at once)
########################################################################################################################
library(pROC)

tiff("GENIE_auc_roc.tiff", width = 1500, height = 1200, units = 'px', res = 300) #create a tiff image
par(pty="s") #remove space between axes and where the plot starts (s=square plotting region)

roc_rf_cv <- pROC::roc(predictor = model_rf$pred$driver, 
                       response = model_rf$pred$obs,
                       levels =c("passenger", "driver"), #set levels for 1st) control vs. 2nd) case
                       direction="<", #set direction < for smaller than or equal median for case vs. control
                       plot=TRUE, 
                       col="green", 
                       lwd=1, #line width of the ROC curve
                       legacy.axes=TRUE, #by default, roc function gives sensitivity vs. specificity, but we need 1-specificity instead
                       grid=TRUE, #add a grid to the plot
                       xlab="False positive rate",
                       ylab="True positive rate",
                       font.lab=2, #font (bold the x-axis and y-axis)
                       font=1, #make tick marks bolder 
                       las=1, #rotate y-axis labels to be 1=horizontal
)

roc_knn_cv <- pROC::roc(predictor = model_knn$pred$driver, 
                        response = model_knn$pred$obs,
                        levels =c("passenger", "driver"), 
                        direction="<", 
                        plot=TRUE, 
                        col="orange", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        add=TRUE)

roc_egb_cv <- pROC::roc(predictor = model_egb$pred$driver, 
                        response = model_egb$pred$obs,
                        levels =c("passenger", "driver"), 
                        direction="<", 
                        plot=TRUE, 
                        col="blueviolet", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        add=TRUE)

roc_svm_cv <- pROC::roc(predictor = model_svm$pred$driver, 
                        response = model_svm$pred$obs,
                        levels =c("passenger", "driver"), 
                        direction="<", 
                        plot=TRUE, 
                        col="chocolate4", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        add=TRUE)

roc_dt_cv <- pROC::roc(predictor = model_dt$pred$driver, 
                       response = model_dt$pred$obs,
                       levels =c("passenger", "driver"), 
                       direction="<",
                       plot=TRUE, 
                       col="blue", 
                       lwd=1, 
                       legacy.axes=TRUE, 
                       add=TRUE)

roc_lr_cv <- pROC::roc(predictor = model_lr$pred$driver, 
                       response = model_lr$pred$obs,
                       levels =c("passenger", "driver"), 
                       direction="<", 
                       plot=TRUE, 
                       col="red", 
                       lwd=1, 
                       legacy.axes=TRUE, 
                       add=TRUE)

roc_mlp_cv <- pROC::roc(predictor = model_mlp$pred$driver, 
                        response = model_mlp$pred$obs,
                        levels =c("passenger", "driver"), 
                        direction="<",
                        plot=TRUE, 
                        col="pink", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        smooth=TRUE,
                        add=TRUE)

abline(a=1, b=-1, col = "black", lwd=0.6)


legend("bottomright",
       legend=c("Random forest (AUC-ROC=0.819)", #add values of AUC-ROC to the legend
                "KNN (AUC-ROC=0.795)",
                "EGB (AUC-ROC=0.783)",
                "SVM (AUC-ROC=0.777)",
                "Decision tree (AUC-ROC=0.757)",
                "Logistic regression (AUC-ROC=0.735)",
                "MLP (AUC-ROC=0.689)",
                "Random chance (AUC-ROC=0.5)"),
       col=c("green", "orange", "blueviolet", "chocolate4", "blue", "red", "pink", "black"), #add colours for each model
       lwd=2, #line width e.g. 10% of normal width 
       pt.lwd = 1, #line width for the points
       cex=0.42, #scale for legend size inside the figure
       pt.cex=0.5, #expansion factor for the points
       text.font=2, #font size of text
       text.width = 0.55, #width of the text (left-right)
       bg="grey" #set colour of the legend background
)

dev.off()

#extract AUC-ROC values of the TRAINING ON GENIE data for all models
roc_lr_cv$auc 
roc_dt_cv$auc 
roc_knn_cv$auc 
roc_egb_cv$auc 
roc_rf_cv$auc 
roc_svm_cv$auc 
roc_mlp_cv$auc


#########################################################################################################################
### EXTRACT AND PLOT AUC-ROC for TESTING ON TCGA data
#########################################################################################################################
library(pROC)

tiff("TCGA_auc_roc.tiff", width = 1500, height = 1200, units = 'px', res = 300) #create a tiff image
par(pty="s") #remove space between axes and where the plot starts (s=square plotting region)

roc_obj_lr <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, lr_predicted_prob$driver,
                        smooth=TRUE, #binormal method of smoothing
                        plot=TRUE, 
                        col="red", 
                        lwd=1, #line width of the ROC curve
                        legacy.axes=TRUE, #by default, roc function gives sensitivity vs. specificity, but we need 1-specificity instead
                        font=1, #make tick marks bolder 
                        las=1, #rotate y-axis labels to be 1=horizontal
                        grid=TRUE, #add a grid to the plot
                        xlab="False positive rate",
                        ylab="True positive rate",
                        font.lab=2, #font (bold the x-axis and y-axis)
)

roc_obj_mlp <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, mlp_predicted_prob$driver,
                         smooth=TRUE, 
                         plot=TRUE, 
                         col="pink", 
                         lwd=1, 
                         legacy.axes=TRUE, 
                         add=TRUE)

roc_obj_dt <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, dt_predicted_prob$driver,
                        smooth=TRUE, 
                        plot=TRUE, 
                        col="blue", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        add=TRUE)

roc_obj_knn <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, knn_predicted_prob$driver, 
                         smooth=TRUE, 
                         plot=TRUE, 
                         col="orange", 
                         lwd=1, 
                         legacy.axes=TRUE, 
                         add=TRUE)

roc_obj_egb <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, egb_predicted_prob$driver,
                         smooth=TRUE, 
                         plot=TRUE, 
                         col="blueviolet", 
                         lwd=1, 
                         legacy.axes=TRUE, 
                         add=TRUE)

roc_obj_rf <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, rf_predicted_prob$driver,
                        smooth=TRUE, 
                        plot=TRUE, 
                        col="green", 
                        lwd=1, 
                        legacy.axes=TRUE, 
                        add=TRUE)

roc_obj_svm <- pROC::roc(under_samp_reduced_TCGA$numeric_labels, svm_predicted_prob$driver,
                         smooth=TRUE, 
                         plot=TRUE, 
                         col="chocolate4", 
                         lwd=1, 
                         legacy.axes=TRUE, 
                         add=TRUE)

abline(a=1, b=-1, col = "black", lwd=0.3)


legend("bottomright",
       legend=c("Logistic regression (AUC-ROC=0.885)", #add values of AUC-ROC to the legend
                "MLP (AUC-ROC=0.876)",
                "Decision tree (AUC-ROC=0.847)",
                "KNN (AUC-ROC=0.833)",
                "EGB (AUC-ROC=0.789)", 
                "Random forest (AUC-ROC=0.785)",
                "SVM (AUC-ROC=0.741)",
                "Random chance (AUC-ROC=0.5)"),
       col=c("red", "pink", "blue", "orange", "blueviolet", "green", "chocolate4", "black"), #add colours for each model
       lwd=2, #line width e.f. 10% of normal width 
       pt.lwd = 1,
       cex=0.42, #scale for legend size inside the figure
       pt.cex=0.5,
       text.font=2, #font size of text
       text.width = 0.55, #width of the text (left-right)
       bg="grey" #set colour of the legend background
)

dev.off()

#extract AUC-ROC values of TRAINING ON GENIE data for all models (smooth=TRUE means with no smoothing vs. smooth=FALSE means with the smoothing function applied)
roc_obj_lr$auc
roc_obj_dt$auc
roc_obj_knn$auc
roc_obj_egb$auc
roc_obj_rf$auc
roc_obj_svm$auc
roc_obj_mlp$auc 




#########################################################################################################################
### PLOTTING A BARCHART WITH PERFORMANCE COMPARISON ON THE BENCHMARK DATASET (MutaGene)
#########################################################################################################################

#generate the table with the performance values
performance_benchmark <- data.frame(Classifier=rep(c("DRIVE", "CHASMplus"),each=5),
                                    Measure=rep(c("Accuracy", "Precision", "Recall", "F1-score", "AUC-ROC"),2),
                                    Score=c(0.852, 0.843, 0.866, 0.854, 0.919, 0.743, 0.672, 0.952, 0.788, 0.872))
#order the rows and column (i.e. as factors)
performance_benchmark$Classifier <- factor(performance_benchmark$Classifier, levels=c("DRIVE", "CHASMplus"))
performance_benchmark$Measure <- factor(performance_benchmark$Measure, levels=c("Accuracy", "Precision", "Recall", "F1-score", "AUC-ROC"))

#generate the barplot using ggplot2
library(ggplot2)

tiff("BENCHMARK_comparison.tiff", width = 1500, height = 1200, units = 'px', res = 300)

bold.text <- element_text(face="bold") #define a bold text for adding to the x and y axes
benchmark_plot <- ggplot(data=performance_benchmark, #data to be used 
                         aes(x=Classifier, y=Score, fill=Measure)) + #define axes and labels
  geom_bar(stat="identity", color="black", position=position_dodge()) + #define outer black lines for each separate bar
  theme(axis.title = bold.text) #make axes bold
benchmark_plot + scale_fill_manual(values=c("chartreuse1","red", "grey", "orange", "slateblue2")) #choose colours for each bar (i.e. performance measure) 

dev.off()
