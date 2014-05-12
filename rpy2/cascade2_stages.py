# -*- coding: utf-8 -*-
"""
Classifier comparisons for 2 stage cascade vs. one-shot all benign all malignant

Created on Mon May 05 15:45:41 2014

@author: Cristina
"""

from rpy2.robjects.packages import importr
caret = importr("caret")
rf = importr("randomForest")
mass = importr("MASS")
pROC = importr("pROC")

utils = importr('utils')
base = importr('base')
stats = importr('stats')
graphics = importr('graphics')
import rpy2.robjects.lib.ggplot2 as ggplot2
from rpy2.robjects.vectors import DataFrame
import rpy2.robjects.packages as rpacks

import rpy2.robjects as robjects


R=10
RFmodel1_AUC = base.numeric(R)
RFmodel1_Sen = base.numeric(R)
RFmodel1_Spec = base.numeric(R)

#for (i in 1:R) 
allRset = utils.read_table("allRset_allFeatures.txt", header=True)
print(allRset.colnames)
print(base.nrow(allRset))

features = utils.read_table("4class_allFeatures.txt", header=True)
print(features.colnames)
print(base.nrow(features))



outcomesetDi  = allRset.rx2("BenignNMaligNAnt")
## the outcome data are needed
partitionsetDi = caret.createDataPartition(y = outcomesetDi, ## the outcome data are needed
                                           p = 0.75,  ## The percentage of data in the training set
                                           list = False) ## The format of the results. 

vecparti = partitionsetDi.rx(True, 1)

nTrainsetDi = base.nrow( allRset.rx(vecparti, True) )
nTestsetDi = robjects.IntVector(range(1))
nTestsetDi.rx[1] = 59
                                   
sTrainsetDi = base.sample(base.nrow(allRset), size=nTrainsetDi ,replace=True)
print(sTrainsetDi)

sTestsetDi = base.sample(base.nrow(allRset), size=nTestsetDi, replace=True)
print(sTestsetDi)

import rpy2.rlike.container as rlc
od = rlc.OrdDict([('BenignNMaligNAnt', allRset.rx2("BenignNMaligNAnt")),
                      ('circularity', allRset.rx2("circularity"))])
                      
dataf = robjects.DataFrame(od)

Set1trainingSel_boot = base.nrow(dataf.rx(sTrainsetDi, True) )
Set1testingSel_boot = base.nrow(dataf.rx(sTestsetDi, True) )


########### train RF
bootControl = caret.trainControl(method = "boot", 
                                number = 10, # number of boostrap iterations
                                savePredictions = True,
                                p = 0.75,
                                classProbs = True,
                                returnResamp = "all",
                                verbose = False,
                                summaryFunction = caret.twoClassSummary)

set1_RFfit = caret.train(BenignNMaligNAnt ~ ., data = Set1trainingSel_boot[,1:10],
                        method = "rf",
                        trControl = bootControl,
                        tuneGrid = robjects.IntVector(range(1,9)),
                        returnData = TRUE,
                        fitBest = TRUE,
                        verbose=FALSE,
                        metric="ROC")
                           

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

RFclassifier = """
compareRF_cascadeTest <- function(R, features, verbose=TRUE) {
  # initialize
  library(caret)
  library(randomForest)
  library(MASS)
  library(mlbench)
  library(pROC)
  
  # R is the number of repetitions (set R = 500 #number of random simulations)
  RFmodel1_AUC = numeric(R)
  RFmodel1_Sen = numeric(R)
  RFmodel1_Spec = numeric(R)
  RFmodel2_AUC = numeric(R)
  RFmodel2_Sen = numeric(R)
  RFmodel2_Spec = numeric(R) 
  RFmodel3_AUC = numeric(R)
  RFmodel3_Sen = numeric(R)
  RFmodel3_Spec = numeric(R) 
  
  cascadeRFmodels_AUC = numeric(R)
  cascadeRFmodels_Sen = numeric(R)
  cascadeRFmodels_Spec = numeric(R)
  oneshotRFmodel4_AUC = numeric(R)
  oneshotRFmodel4_Sen = numeric(R)
  oneshotRFmodel4_Spec = numeric(R)
  pvaltestROC_RFcascadevsRFoneshot = numeric(R)
  deltaAUC_RFcascadevsRFoneshot = numeric(R)

  # create data.frame to hold resamples of cascade
  cascadeCollected <- 1:10
  names(cascadeCollected) = c("caseID", "labels", "pred1", "obs1", "P1", "N1", "pred2",    "obs2",    "P2", "N2")
  
  # start R iterations
  for (i in 1:R) {
    #####################################
    # 1) Set bootstrap resamples sizes
    #####################################
    set <-  na.omit(features[c(10:87)])
    setDi <- set[c(1)]
    outcomesetDi  <- setDi$BenignNMaligNAnt
    partitionsetDi <- createDataPartition(y = outcomesetDi, ## the outcome data are needed
                                         p = .75, ## The percentage of data in the training set
                                         list = FALSE) ## The format of the results. 
    
    nTrainsetDi <-  length( setDi[ partitionsetDi ,] )
    nTestsetDi <-  length( setDi[-partitionsetDi ,] )
    
    ## 2) sampled with replacement from the complete cross validation data set 5000 times
    sTrainsetDi <- sample(nrow(set),size=nTrainsetDi ,replace=T)
    sTestsetDi<- sample(nrow(set),size=nTestsetDi, replace=T)
    
    ################################################
    # 2) Select subsets of features correspondingly
    ################################################
    set1_selfeatures <- set[c("BenignNMaligNAnt",
                                      "circularity",
                                      "Tpeak.inside",
                                      "Slope_ini.inside",
                                      "SER.inside",
                                      "SER.countor",
                                      "Kpeak.inside",
                                      "iiMin_change_Variance_uptake",
                                      "iMax_Variance_uptake",
                                      "Kpeak.countor")]
    # generate case IDS to recollect after
    set1_selfeatures$caseID = rownames(set1_selfeatures)
    set1_selfeatures$labels <- set1_selfeatures$BenignNMaligNAnt
    
    Set1trainingSel.boot <- set1_selfeatures[sTrainsetDi ,]
    Set1testingSel.boot <- set1_selfeatures[sTestsetDi ,]
        
    M<-subset(Set1trainingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="massM")
    ifelse( M$BenignNMaligNAnt == "massB", "mass", "mass") -> M$BenignNMaligNAnt
    N<-subset(Set1trainingSel.boot, BenignNMaligNAnt=="nonmassB" | BenignNMaligNAnt=="nonmassM")
    ifelse( N$BenignNMaligNAnt == "nonmassB", "nonmass", "nonmass") -> N$BenignNMaligNAnt
    
    Set1trainingSel.boot = data.frame(rbind(M,N))
    
    M<-subset(Set1testingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="massM")
    ifelse( M$BenignNMaligNAnt == "massB", "mass", "mass") -> M$BenignNMaligNAnt
    N<-subset(Set1testingSel.boot, BenignNMaligNAnt=="nonmassB" | BenignNMaligNAnt=="nonmassM")
    ifelse( N$BenignNMaligNAnt == "nonmassB", "nonmass", "nonmass") -> N$BenignNMaligNAnt
    
    Set1testingSel.boot <- data.frame(rbind(M,N))

    ################## 
    ### Select subsets of features correspondingly
    set2_selfeatures <- set[c("BenignNMaligNAnt","Kpeak.countor",
                              "beta.countor",
                              "irregularity",
                              "washoutRate.contour",
                              "UptakeRate.contour",
                              "SER.countor",
                              "Tpeak.countor",
                              "A.countor",
                              "min_F_r_i",
                              "maxCr.contour",
                              "circularity",
                              "beta.inside",
                              "UptakeRate.inside",
                              "texture_energy_threeQuaRad",
                              "iAUC1.countor",
                              "texture_ASM_threeQuaRad",
                              "Kpeak.inside",
                              "alpha.inside",
                              "Vr_increasingRate.inside",
                              "maxVr.inside",
                              "texture_ASM_quarterRad",
                              "texture_energy_quarterRad",
                              "maxCr.inside",
                              "max_F_r_i",
                              "A.inside",
                              "texture_energy_zero",
                              "iAUC1.inside",
                              "texture_contrast_quarterRad",
                              "Tpeak.inside",
                              "peakCr.contour",
                              "texture_contrast_halfRad",
                              "Vr_post_1.inside",
                              "alpha.countor",
                              "texture_homogeneity_quarterRad",
                              "texture_contrast_threeQuaRad",
                              "Slope_ini.countor",
                              "texture_ASM_zero",
                              "texture_energy_halfRad",
                              "texture_ASM_halfRad",
                              "texture_contrast_zero",
                              "Vr_decreasingRate.contour",
                              "SER.inside",
                              "mean_F_r_i",
                              "Vr_increasingRate.contour",
                              "var_F_r_i")]
    
    # generate case IDS to recollect after
    set2_selfeatures$caseID = rownames(set2_selfeatures)
    set2_selfeatures$labels <- set2_selfeatures$BenignNMaligNAnt
    
    Set2trainingSel.boot <- set2_selfeatures[sTrainsetDi ,]
    Set2testingSel.boot <- set2_selfeatures[sTestsetDi ,]
    
    M<-subset(Set2trainingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="massM")
    ifelse( M$BenignNMaligNAnt == "massB", "NC", "C") -> M$BenignNMaligNAnt    
    Set2trainingSel.boot = data.frame(M)
    
    M<-subset(Set2testingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="massM")
    ifelse( M$BenignNMaligNAnt == "massB", "NC", "C") -> M$BenignNMaligNAnt    
    Set2testingSel.boot = data.frame(M)
    
    ################## 
    ### Select subsets of features correspondingly
    set3_selfeatures <- set[c("BenignNMaligNAnt",
                              "Slope_ini.countor",
                              "edge_sharp_std",
                              "UptakeRate.contour",
                              "var_F_r_i",
                              "maxCr.contour",
                              "Tpeak.countor",
                              "alpha.countor",
                              "max_F_r_i",
                              "iiMin_change_Variance_uptake",
                              "max_RGH_mean",
                              "texture_correlation_quarterRad",
                              "max_RGH_var",
                              "texture_contrast_threeQuaRad",
                              "max_RGH_var_k",
                              "washoutRate.contour")]
    
    # generate case IDS to recollect after
    set3_selfeatures$caseID = rownames(set3_selfeatures)
    set3_selfeatures$labels <- set3_selfeatures$BenignNMaligNAnt
    
    Set3trainingSel.boot <- set3_selfeatures[sTrainsetDi ,]
    Set3testingSel.boot <- set3_selfeatures[sTestsetDi ,]
    
    N<-subset(Set3trainingSel.boot, BenignNMaligNAnt=="nonmassB" | BenignNMaligNAnt=="nonmassM")
    ifelse( N$BenignNMaligNAnt == "nonmassB", "NC", "C") -> N$BenignNMaligNAnt    
    Set3trainingSel.boot = data.frame(N)
    
    N<-subset(Set3testingSel.boot, BenignNMaligNAnt=="nonmassB" | BenignNMaligNAnt=="nonmassM")
    ifelse( N$BenignNMaligNAnt == "nonmassB", "NC", "C") -> N$BenignNMaligNAnt    
    Set3testingSel.boot = data.frame(N)
    
    #####################################
    # 3) Train Random Forest Classifiers
    #####################################
    bootControl <- trainControl(method = "boot", 
                                number = 10, # number of boostrap iterations
                                savePredictions = TRUE,
                                p = 0.75,
                                classProbs = TRUE,
                                returnResamp = "all",
                                verbose = FALSE,
                                summaryFunction = twoClassSummary)
   
    ########## RF set1
    RFGrid <- expand.grid(.mtry=c(1:9) )
    #set.seed(122)
    set1_RFfit <- train(as.factor(BenignNMaligNAnt) ~ ., data = Set1trainingSel.boot[,1:10],
                        method = "rf",
                        trControl = bootControl,
                        tuneGrid = RFGrid,
                        returnData = TRUE,
                        fitBest = TRUE,
                        verbose=FALSE,
                        metric="ROC")
    if(verbose) print(set1_RFfit$finalModel)
    # Classify Training data to get training error and accuracy
    RF1_train <- confusionMatrix(predict(set1_RFfit , newdata = Set1trainingSel.boot[,1:10]), Set1trainingSel.boot$BenignNMaligNAnt )
    RF1_test <- confusionMatrix(predict(set1_RFfit , newdata = Set1testingSel.boot[,1:10]), Set1testingSel.boot$BenignNMaligNAnt )
  if(verbose) print(RF1_test)
    
    ## compute classifier stage 1 performance    
    RFmodel1 <- list(RFmodel1 = set1_RFfit)
    RFmodel1_predValues <- extractPrediction(RFmodel1, testX= Set1trainingSel.boot[,2:10], testY= Set1trainingSel.boot[,1])
    RFmodel1_predValues <- extractPrediction(RFmodel1, testX= Set1testingSel.boot[,2:10], testY= Set1testingSel.boot[,1])
    if(verbose) table(RFmodel1_predValues)
    RFmodel1_probValues <- extractProb(RFmodel1, testX=Set1testingSel.boot[,2:10], testY=Set1testingSel.boot[,1])
    RFmodel1_testProbs <- subset(RFmodel1_probValues, dataType=="Test")
    
    # Edit probabilities so that can do pooled ROC
    RFmodels_cascade <- RFmodel1_probValues
    RFmodels_cascade$labels <- RFmodels_cascade$obs
    ifelse( RFmodels_cascade$obs == "mass", "P", "N") -> RFmodels_cascade$obs
    ifelse( RFmodels_cascade$pred == "mass", "P", "N") -> RFmodels_cascade$pred
    colnames(RFmodels_cascade)[1] <- "P"
    colnames(RFmodels_cascade)[2] <- "N"
    
    ## plot
    RFmodel1_ROC <- roc(RFmodel1_testProbs$obs, RFmodel1_testProbs$mass,
                             main="Classifier for mass vs. nonmass",
                             percent=TRUE,
                             ci = TRUE,
                             of = "se", 
                             smooth.method="binormal",
                             sp = seq(0, 100, 10))
    if(verbose) print(RFmodel1_ROC) 
  
    ########### RF set2
    RFGrid <- expand.grid(.mtry=c(1:10,15,20,25,30,35,40) )
    #set.seed(122)
    set2_RFfit <- train(as.factor(BenignNMaligNAnt) ~ ., data = Set2trainingSel.boot[1:46],
                        method = "rf",
                        trControl = bootControl,
                        tuneGrid = RFGrid,
                        returnData = TRUE,
                        fitBest = TRUE,
                        verbose=FALSE,
                        metric="ROC")
    if(verbose) print(set2_RFfit$finalModel)
    
    # Classify Training data to get training error and accuracy
    RF2_train <- confusionMatrix(predict(set2_RFfit , newdata = Set2trainingSel.boot), Set2trainingSel.boot$BenignNMaligNAnt )
    RF2_test <- confusionMatrix(predict(set2_RFfit , newdata = Set2testingSel.boot), Set2testingSel.boot$BenignNMaligNAnt )
    if(verbose) print(RF2_test)
    
    ###### Compute ROC
    RFmodel2 <- list(RFmodel2 = set2_RFfit)
    RFmodel2_predValues <- extractPrediction(RFmodel2, testX= Set2trainingSel.boot[,2:46], testY= Set2trainingSel.boot[,1])
    RFmodel2_predValues <- extractPrediction(RFmodel2, testX= Set2testingSel.boot[,2:46], testY= Set2testingSel.boot[,1])
    if(verbose) table(RFmodel2_predValues)
    
    RFmodel2_probValues <- extractProb(RFmodel2, testX=Set2testingSel.boot[,2:46], testY=Set2testingSel.boot[,1])
    RFmodel2_testProbs <- subset(RFmodel2_probValues, dataType=="Test")
    
    # Edit probabilities so that can do pooled ROC
    ifelse( RFmodel2_probValues$obs == "C", "P", "N") -> RFmodel2_probValues$obs
    ifelse( RFmodel2_probValues$pred == "C", "P", "N") -> RFmodel2_probValues$pred
    colnames(RFmodel2_probValues)[1] <- "P"
    colnames(RFmodel2_probValues)[2] <- "N"
    
    RFmodel2_ROC <- roc(RFmodel2_testProbs$obs, RFmodel2_testProbs$C,
                             percent=TRUE, 
                             ci = TRUE,
                             of = "se", 
                             smooth.method="binormal",
                             sp = seq(0, 100, 10))
    if(verbose) print(RFmodel2_ROC)
         
    ########### RF set3
    RFGrid <- expand.grid(.mtry=c(1:10) )
    #set.seed(122)
    set3_RFfit <- train(as.factor(BenignNMaligNAnt) ~ ., data = Set3trainingSel.boot[1:16],
                        method = "rf",
                        trControl = bootControl,
                        tuneGrid = RFGrid,
                        returnData = TRUE,
                        fitBest = TRUE,
                        verbose=FALSE,
                        metric="ROC")
    if(verbose)  print(set3_RFfit$finalModel)
    
    # Classify Training data to get training error and accuracy
    RF3_train <- confusionMatrix(predict(set3_RFfit , newdata = Set3trainingSel.boot), Set3trainingSel.boot$BenignNMaligNAnt )
    RF3_test <- confusionMatrix(predict(set3_RFfit , newdata = Set3testingSel.boot), Set3testingSel.boot$BenignNMaligNAnt )
    if(verbose) print(RF3_test)
    
    ###### Compute ROC
    RFmodel3 <- list(RFmodel3 = set3_RFfit)
    RFmodel3_predValues <- extractPrediction(RFmodel3, testX= Set3trainingSel.boot[,2:16], testY= Set3trainingSel.boot[,1])
    RFmodel3_predValues <- extractPrediction(RFmodel3, testX= Set3testingSel.boot[,2:16], testY= Set3testingSel.boot[,1])
    if(verbose) table(RFmodel3_predValues)
    
    RFmodel3_probValues <- extractProb(RFmodel3, testX=Set3testingSel.boot[,2:16], testY=Set3testingSel.boot[,1])
    RFmodel3_testProbs <- subset(RFmodel3_probValues, dataType=="Test")
    
    # Edit probabilities so that can do pooled ROC
    ifelse( RFmodel3_probValues$obs == "C", "P", "N") -> RFmodel3_probValues$obs
    ifelse( RFmodel3_probValues$pred == "C", "P", "N") -> RFmodel3_probValues$pred
    colnames(RFmodel3_probValues)[1] <- "P"
    colnames(RFmodel3_probValues)[2] <- "N"
    
    RFmodel3_ROC <- roc(RFmodel3_testProbs$obs, RFmodel3_testProbs$C,
                             percent=TRUE, 
                             ci = TRUE,
                             of = "se", 
                             smooth.method="binormal",
                             sp = seq(0, 100, 10))
    if(verbose) {
      print(RFmodel3_ROC)
      plot.roc(RFmodel1_ROC,col="#008600")
      par(new=TRUE)
      plot.roc(RFmodel2_ROC, col="#860000")
      par(new=TRUE)
      plot.roc(RFmodel3_ROC, col="#000086")
      legend("bottomright", legend=c("mass/nonmass", "BM only mass", "BM only nonmass"), col=c("#008600","#860000", "#000086"), lwd=2)
    } 
      
    #################
    # 4) Results
    #################
    ### Print results
    if(verbose) {
      print(i)
      print("2class_MassorNonmass:")
      print(RFmodel1_ROC$auc)
    }
    RFmodel1_AUC[i] = RFmodel1_ROC$auc
    
    ####   postResample Sensitivity & Specificity  
    Set1_PredictResTest <- predict(set1_RFfit, newdata = Set1testingSel.boot)
    RFmodel1_Sen[i] = sensitivity(Set1_PredictResTest, as.factor(Set1testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel1_Sen[i])
    RFmodel1_Spec[i] = specificity(Set1_PredictResTest, as.factor(Set1testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel1_Spec[i])
    
    ### Print results
    if(verbose) {
      print("2classes_mass:")
      print(RFmodel2_ROC$auc)
    }
    RFmodel2_AUC[i] = RFmodel2_ROC$auc
    
    ####   postResample Sensitivity & Specificity  
    Set2_PredictResTest <- predict(set2_RFfit, newdata = Set2testingSel.boot)
    RFmodel2_Sen[i] = sensitivity(Set2_PredictResTest, as.factor(Set2testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel2_Sen[i])
    RFmodel2_Spec[i] = specificity(Set2_PredictResTest, as.factor(Set2testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel2_Spec[i])
    
    ### Print results
    if(verbose) {
      print("2classes_nonmass:")
      print(RFmodel3_ROC$auc)
    }  
    RFmodel3_AUC[i] = RFmodel3_ROC$auc
    
    ####   postResample Sensitivity & Specificity	
    Set3_PredictResTest <- predict(set3_RFfit, newdata = Set3testingSel.boot)
    RFmodel3_Sen[i] = sensitivity(Set3_PredictResTest, as.factor(Set3testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel3_Sen[i])
    RFmodel3_Spec[i] = specificity(Set3_PredictResTest, as.factor(Set3testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(RFmodel3_Spec[i])
    
    #########################################################
    # 5) train BorM binary and compare to other classifier
    #########################################################
    ## 1) Set resamples sizes
    set4_selfeatures <- set[c("BenignNMaligNAnt",
                              "circularity",
                              "Tpeak.inside",
                              "Slope_ini.inside",
                              "SER.inside",
                              "SER.countor",
                              "Kpeak.inside",
                              "iiMin_change_Variance_uptake",
                              "iMax_Variance_uptake",
                              "Kpeak.countor")]
    
    Set4trainingSel.boot <- set4_selfeatures[sTrainsetDi ,]
    Set4testingSel.boot <- set4_selfeatures[sTestsetDi ,]
    
    C<-subset(Set4trainingSel.boot, BenignNMaligNAnt=="massM" | BenignNMaligNAnt=="nonmassM")
    ifelse( C$BenignNMaligNAnt == "massM", "C", "C") -> C$BenignNMaligNAnt
    NC<-subset(Set4trainingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="nonmassB")
    ifelse( NC$BenignNMaligNAnt == "massB", "NC", "NC") -> NC$BenignNMaligNAnt
    
    Set4trainingSel.boot = data.frame(rbind(C,NC))
    
    C<-subset(Set4testingSel.boot, BenignNMaligNAnt=="massM" | BenignNMaligNAnt=="nonmassM")
    ifelse( C$BenignNMaligNAnt == "massM", "C", "C") -> C$BenignNMaligNAnt
    NC<-subset(Set4testingSel.boot, BenignNMaligNAnt=="massB" | BenignNMaligNAnt=="nonmassB")
    ifelse( NC$BenignNMaligNAnt == "massB", "NC", "NC") -> NC$BenignNMaligNAnt
    
    Set4testingSel.boot <- data.frame(rbind(C,NC))
    
    ##################
    ## 3) For each resampling two AUC ROC curves are constructed for the two feature sets to be compared.
    bootControl <- trainControl(method = "boot", 
                                number = 10, # number of boostrap iterations
                                savePredictions = TRUE,
                                p = 0.75,
                                classProbs = TRUE,
                                returnResamp = "all",
                                verbose = FALSE,
                                summaryFunction = twoClassSummary)
    #OUTPUT
    ########## RF 1
    RFGrid <- expand.grid(.mtry=c(1:8) )
    #set.seed(1)
    set4_RFfit <- train(as.factor(BenignNMaligNAnt) ~ ., data = Set4trainingSel.boot,
                        method = "rf",
                        trControl = bootControl,
                        tuneGrid = RFGrid,
                        returnData = TRUE,
                        fitBest = TRUE,
                        verbose=FALSE,
                        metric="ROC")
    if(verbose) print(set4_RFfit$finalModel)
    
    # Classify Training data to get training error and accuracy
    massRF4_test <- confusionMatrix(predict(set4_RFfit , newdata = Set4testingSel.boot), Set4testingSel.boot$BenignNMaligNAnt )
    if(verbose) print(massRF4_test)    
    
    #############################################
    # 6) Final Results: pooling of data for ROC
    #############################################
    # combining RF probs (initial experiment)
    # RFmodels_probValues = rbind(RFmodel1_probValues,RFmodel2_probValues,RFmodel3_probValues)
    # RFmodels_testProbs <- subset(RFmodels_probValues, dataType=="Test")
    # Now use casetest_cascade to fed them to 2nd stage classifiers based on predictions of 1st classifier
    # Record their new predictions
    ##################
    RFmodel1 <- list(RFmodel1 = set1_RFfit)
    
    k=1
    probValues <- extractProb(RFmodel1, testX=Set1testingSel.boot[k,2:10], testY=Set1testingSel.boot[k,1])
    k_probValues <- subset(probValues, dataType=="Test")
    k_caseID = Set1testingSel.boot[k,11]
    casestest = data.frame(Set1testingSel.boot[Set1testingSel.boot$caseID == k_caseID,11:12])
    casestest$pred1 = k_probValues$pred
    casestest$obs1 = k_probValues$obs
    casestest$P1 = k_probValues$mass
    casestest$N1 = k_probValues$nonmass
    
    # initialize dataframe
    casestest_cascade = casestest
    
    for (k in 2:nrow(Set1testingSel.boot)) {
      probValues <- extractProb(RFmodel1, testX=Set1testingSel.boot[k,2:10], testY=Set1testingSel.boot[k,1])
      k_probValues <- subset(probValues, dataType=="Test")
      
      k_caseID = Set1testingSel.boot[k,11]
      casestest = data.frame(Set1testingSel.boot[Set1testingSel.boot$caseID == k_caseID,11:12])
      casestest$pred1 = k_probValues$pred
      casestest$obs1 = k_probValues$obs
      casestest$P1 = k_probValues$mass
      casestest$N1 = k_probValues$nonmass
      #append to dataframe
      casestest_cascade = rbind(casestest_cascade, casestest)
    }
    
    # instantiate the 2nd stage RF classifiers
    RFmodelM <- list(RFmodelM = set2_RFfit)
    RFmodelNonM <- list(RFmodelNonM = set3_RFfit)
    
    #subseting to fed into 2nd stage
    M.boot <- subset(casestest_cascade, pred1=="mass")
    NonM.boot<- subset(casestest_cascade, pred1=="nonmass")
    
    k=1
    M_caseID = M.boot[k,1]
    caseM.boot = data.frame(set2_selfeatures[set2_selfeatures$caseID == M_caseID,])
    ifelse( caseM.boot$BenignNMaligNAnt == "massB", "NC", "C") -> caseM.boot$BenignNMaligNAnt    
    SetMtestingSel.boot = data.frame(caseM.boot)
      
    probValuesM <- extractProb(RFmodelM, testX=SetMtestingSel.boot[,2:46], testY=SetMtestingSel.boot[,1])
    k_probValuesM <- subset(probValuesM, dataType=="Test")
    
    M.boot_frame = data.frame(M.boot[k,])
    M.boot_frame$pred2 = k_probValuesM$pred
    M.boot_frame$obs2 = k_probValuesM$obs
    M.boot_frame$P2 = k_probValuesM$C
    M.boot_frame$N2 = k_probValuesM$NC
    
    # initialize dataframe
    M.boot_cascade = M.boot_frame
    SetMtestingSel.boot_cascade = data.frame(caseM.boot)
    
    for (k in 2:nrow(M.boot)) {
      M_caseID = M.boot[k,1]
      caseM.boot = data.frame(set2_selfeatures[set2_selfeatures$caseID == M_caseID,])
      ifelse( caseM.boot$BenignNMaligNAnt == "massB", "NC", "C") -> caseM.boot$BenignNMaligNAnt    
      SetMtestingSel.boot = data.frame(caseM.boot)
      
      probValuesM <- extractProb(RFmodelM, testX=SetMtestingSel.boot[,2:46], testY=SetMtestingSel.boot[,1])
      k_probValuesM <- subset(probValuesM, dataType=="Test")
      
      M.boot_frame = data.frame(M.boot[k,])
      M.boot_frame$pred2 = k_probValuesM$pred 
      M.boot_frame$obs2 = k_probValuesM$obs
      M.boot_frame$P2 = k_probValuesM$C
      M.boot_frame$N2 = k_probValuesM$NC
      
      #append to dataframe
      M.boot_cascade = rbind(M.boot_cascade, M.boot_frame)
      SetMtestingSel.boot_cascade = rbind(SetMtestingSel.boot_cascade, caseM.boot)
    }
    
    # Do non-masses frame
    k=1
    NonM_caseID = NonM.boot[k,1]
    caseNonM.boot = data.frame(set3_selfeatures[set3_selfeatures$caseID == NonM_caseID,])
    ifelse( caseNonM.boot$BenignNMaligNAnt == "massB", "NC", "C") -> caseNonM.boot$BenignNMaligNAnt    
    SetNonMtestingSel.boot = data.frame(caseNonM.boot)
    
    probValuesNonM <- extractProb(RFmodelNonM, testX=SetNonMtestingSel.boot[,2:16], testY=SetNonMtestingSel.boot[,1])
    k_probValuesNonM <- subset(probValuesNonM, dataType=="Test")
    
    NonM.boot_frame = data.frame(NonM.boot[k,])
    NonM.boot_frame$pred2 = k_probValuesNonM$pred
    NonM.boot_frame$obs2 = k_probValuesNonM$obs
    NonM.boot_frame$P2 = k_probValuesNonM$C
    NonM.boot_frame$N2 = k_probValuesNonM$NC
    
    # initialize dataframe
    NonM.boot_cascade = NonM.boot_frame
    
    for (k in 2:nrow(NonM.boot)) {
      NonM_caseID = NonM.boot[k,1]
      caseNonM.boot = data.frame(set3_selfeatures[set3_selfeatures$caseID == NonM_caseID,])
      ifelse( caseNonM.boot$BenignNMaligNAnt == "nonmassB", "NC", "C") -> caseNonM.boot$BenignNMaligNAnt    
      SetNonMtestingSel.boot = data.frame(caseNonM.boot)
      
      probValuesNonM <- extractProb(RFmodelNonM, testX=SetNonMtestingSel.boot[,2:16], testY=SetNonMtestingSel.boot[,1])
      k_probValuesNonM <- subset(probValuesNonM, dataType=="Test")
      
      NonM.boot_frame = data.frame(NonM.boot[k,])
      NonM.boot_frame$pred2 = k_probValuesNonM$pred 
      NonM.boot_frame$obs2 = k_probValuesNonM$obs
      NonM.boot_frame$P2 = k_probValuesNonM$C
      NonM.boot_frame$N2 = k_probValuesNonM$NC
      
      #append to dataframe
      NonM.boot_cascade = rbind(NonM.boot_cascade, NonM.boot_frame)
    }
            
    ##################  
    # generate an empy dataframe to copy ouput probabilities of stage 1 classifier
    RFcascade_stage1 <- as.data.frame(setNames(replicate(10,numeric(nrow(RFmodels_cascade)), simplify = F), colnames(M.boot_cascade)))
    RFcascade_stage1$caseID = 1:nrow(RFmodels_cascade)
    RFcascade_stage1$labels = "mass"
    RFcascade_stage1$pred1 = RFmodels_cascade$pred
    RFcascade_stage1$obs1 = RFmodels_cascade$obs
    RFcascade_stage1$P1 = RFmodels_cascade$P
    RFcascade_stage1$N1 = RFmodels_cascade$N
    RFcascade_stage1$pred2 = RFmodels_cascade$pred
    RFcascade_stage1$obs2 = RFmodels_cascade$obs
    RFcascade_stage1$P2 = RFmodels_cascade$P
    RFcascade_stage1$N2 = RFmodels_cascade$N
    
    # 3) Final Results: pooling of data for ROC
    RFmodels_probValues = rbind(M.boot_cascade,NonM.boot_cascade)
    RFmodels_probValuesPooled = rbind(M.boot_cascade,NonM.boot_cascade, RFcascade_stage1)
        
    RFmodel4 <- list(RFmodel4 = set4_RFfit)
    RFmodel4_predValues <- extractPrediction(RFmodel4, testX= Set4testingSel.boot[,-1:-1], testY= Set4testingSel.boot[,1])
    table(RFmodel4_predValues)
    
    RFmodel4_probValues <- extractProb(RFmodel4, testX=Set4testingSel.boot[,-1:-1], testY=Set4testingSel.boot[,1])
    RFmodel4_testProbs <- subset(RFmodel4_probValues, dataType=="Test")
    
    ########################################## ###### ###### 
    ###### Compute ROC
    RFmodels_ROC <- roc(RFmodels_probValuesPooled$obs2, RFmodels_probValuesPooled$P2,
                             main="Cascade Classifier for mass and B or M (pooled data)",
                             percent=TRUE,
                             ci = TRUE,
                             of = "se", 
                             smooth.method="binormal",
                             sp = seq(0, 100, 10))
    
    if(verbose) print(RFmodels_ROC)

    RFmodel4_ROC <- roc(RFmodel4_testProbs$obs, RFmodel4_testProbs$C,
                             percent=TRUE,
                             ci = TRUE,
                             of = "se", 
                             smooth.method="binormal",
                             sp = seq(0, 100, 10))
    if(verbose) print(RFmodel4_ROC)
    
    ###################### Print results
    if(verbose) {
      print(i)
      print("Cascade Results for cascadeRFmodels_ROC:")
      print(RFmodels_ROC$auc)
    }
    cascadeRFmodels_AUC[i] = RFmodels_ROC$auc
    
    ####   postResample Sensitivity & Specificity 
    SetMass_PredictResTest <- predict(set2_RFfit, newdata = SetMtestingSel.boot_cascade[,1:46])
    cascadeRFmodels_Sen[i] = sensitivity(SetMass_PredictResTest, as.factor(SetMtestingSel.boot_cascade$BenignNMaligNAnt))
    if(verbose) print(cascadeRFmodels_Sen[i])
    cascadeRFmodels_Spec[i] = specificity(SetMass_PredictResTest, as.factor(SetMtestingSel.boot_cascade$BenignNMaligNAnt))
    if(verbose) print(cascadeRFmodels_Spec[i])
    
    if(verbose) {
      print("Results for One-shot:")
      print(RFmodel4_ROC$auc)
      oneshotRFmodel4_AUC[i] = RFmodel4_ROC$auc
    }
    
    ####   postResample Sensitivity & Specificity  
    Set4_PredictResTest <- predict(set4_RFfit, newdata = Set4testingSel.boot)
    oneshotRFmodel4_Sen[i] = sensitivity(Set4_PredictResTest, as.factor(Set4testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(oneshotRFmodel4_Sen[i])
    oneshotRFmodel4_Spec[i] = specificity(Set4_PredictResTest, as.factor(Set4testingSel.boot$BenignNMaligNAnt))
    if(verbose) print(oneshotRFmodel4_Spec[i])
    
    # do t-test for AUC
    pvaltestROC <- roc.test(RFmodels_ROC, RFmodel4_ROC)
    pvaltestROC_RFcascadevsRFoneshot[i] = pvaltestROC$p.value
    
    if(verbose) {
      plot.roc(RFmodels_ROC, col="#860000")
      par(new=TRUE)
      plot.roc(RFmodel4_ROC, col="#008600")
      legend("bottomright", legend=c("cascade MvsN/CvsNC", "BM one-shot"), col=c("#860000", "#008600"), lwd=2)
      text(50, 50, labels=paste("p-value =", format.pval(pvaltestROC$p.value)), adj=c(0, .5))
    }
    
    ###########    
    ## 4) The difference in AUC ROC was computed. Resampling 5000 times resulted in 1000 values for deltaAUC
    deltaAUC_RFcascadevsRFoneshot[i] = cascadeRFmodels_AUC[i] - oneshotRFmodel4_AUC[i]
    
    # collect some resamples to evaluate retrospectively
    if( i %% 2 == 0){
      cascadeCollected = rbind(cascadeCollected, RFmodels_probValues) 
    }
  }
  
  ######## Finally write output  
  output<-list(RFmodel1_AUC=RFmodel1_AUC, RFmodel1_Sen=RFmodel1_Sen, RFmodel1_Spec=RFmodel1_Spec, 
               RFmodel2_AUC=RFmodel2_AUC, RFmodel2_Sen=RFmodel2_Sen, RFmodel2_Spec=RFmodel2_Spec, 
               RFmodel3_AUC=RFmodel3_AUC, RFmodel3_Sen=RFmodel3_Sen, RFmodel3_Spec=RFmodel3_Spec, 
               cascadeRFmodels_AUC=cascadeRFmodels_AUC, cascadeRFmodels_Sen=cascadeRFmodels_Sen, 
               cascadeRFmodels_Spec=cascadeRFmodels_Spec, oneshotRFmodel4_AUC=oneshotRFmodel4_AUC,
               oneshotRFmodel4_Sen=oneshotRFmodel4_Sen, oneshotRFmodel4_Spec=oneshotRFmodel4_Spec,
               pvaltestROC_RFcascadevsRFoneshot=pvaltestROC_RFcascadevsRFoneshot,
               deltaAUC_RFcascadevsRFoneshot=deltaAUC_RFcascadevsRFoneshot,
               cascadeCollected=cascadeCollected)

  return(output)  
}
 """
 
 rfun = SignatureTranslatedAnonymousPackage(RFclassifier, "rfun")
                           