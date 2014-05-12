subset_feature_selection <- function(filename_data, seed) {
  library(caret)
  library(mlbench)
  library(pROC)
  
  # Using RF functions
  # read data into a dataframe called "features "
  features <- read.table(filename_data, header=T)
  
  # Only select 10th thru 93th variables: excluding "N_interNAlvoxels"  "n_VOI_points" 
  allfeatures <- na.omit(features[c(10:87)]) #[c(1,5:10)] # select 1st and 5th thru 10th variables
  
  ### Pre-p[rocess
  predictors <- allfeatures
  outcome <- na.omit(allfeatures[c(1)]) #[c(1,5:10)] # select 1st and 5th thru 10th variables
  outcome  <- outcome$BenignNMaligNAnt
  
  set.seed(2)
  inTrain <- createDataPartition(y = outcome, ## the outcome data are needed
                                 p = .75, ## The percentage of data in the training set
                                 list = FALSE) ## The format of the results. 
  
  training <- predictors[ inTrain,]
  testing <- predictors[-inTrain,]
  nrow( training )
  nrow( testing )
  
  ############ Recursive Feature Selection via caret 
  # fit models with subset sizes of 1:10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70
  subsets <- c(1:10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70)
  
  # create control object for Controlling the Feature Selection Algorithms
  # Right now performs 10 repeated cross validation
  RFctrl <- rfeControl(functions = rfFuncs, 
                       method = "repeatedcv", 
                       repeats = 10,
                       verbose = FALSE,
                       returnResamp = "all")
  
  set.seed(seed)
  # Run recursive feature selection (RFE) algorithm
  rfSelProfile <- rfe(predictors[,2:78], outcome,
                      sizes = subsets, rfeControl = RFctrl)
  
  print(rfSelProfile )
  
  # The predictors function can be used to get a text string of variable names that were picked in      the final model. 
  # The model can be used to get best subset and predictions for future or test samples.
  print(rfSelProfile$bestSubset)
  print(rfSelProfile$optVariables)
  
  # Also the resampling results are stored in the sub-object 
  # and can be used with several lattice functions. Univariate lattice functions (densityplot, histogram) 
  # can be used to plot the resampling distribution while bivariate functions (xyplot, stripplot) 
  # can be used to plot the distributions for different subset sizes.
  #head(rfSelProfile$resample)
  
  # plot to visualize the results. 
  plot(rfSelProfile, type = c("g", "o"))
  
  ## Create dataframe with one selected features
  selfeatures = data.frame(training[,c("BenignNMaligNAnt",rfSelProfile$optVariables)])
  
  ################
  ## For picking subset sizes:
  ## Maximize Accuracy
  performance <- data.frame(Accuracy = rfSelProfile$results$Accuracy,
                            Variables = rfSelProfile$results$Variables)
  
  ## Percent Loss in performance (positive)
  performance$PctLoss <- (max(performance$Accuracy ) - performance$Accuracy )/max(performance$Accuracy )*100
  
  plot(performance$Variables , performance$Accuracy, type="p",  col="blue", xlab="Variables", ylab="Accuracy")
  lines(performance$Variables , performance$Accuracy, col="blue") 
  grid(10, 10, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  Axis(at =  seq(0, 80, 5), side=1, labels = TRUE)
  
  absoluteBest <- pickSizeBest(performance, metric = "Accuracy", maximize = TRUE)
  within5Pct <- pickSizeTolerance(performance, metric = "Accuracy", maximize = TRUE)
  
  cat("numerically optimal:", performance$Accuracy[performance$Variables==absoluteBest ],"Accuracy with subset",absoluteBest, "\n")
  cat("Accepting a 1.5% Accuracy loss:", performance$Accuracy[performance$Variables==within5Pct ],"Accuracy with subset",within5Pct, "\n")
  
  plot.new()
  plot(performance$Variables , performance$PctLoss, type="p",  col="blue", xlab="Variables", ylab="Accuracy % Loss")
  lines(performance$Variables , performance$PctLoss, col="blue") 
  ### Add those points to plot
  points(absoluteBest, performance$PctLoss[performance$Variables==absoluteBest], type="p", col="red", bg="red", pch=22, lwd=1)
  points(within5Pct, performance$PctLoss[performance$Variables==within5Pct], type="p", col="black", bg="black", pch=25, lwd=1)
  # Organize plot
  grid(10, 10, col = "lightgray", lty = "dotted", lwd = par("lwd"))
  Axis(at =  seq(0, 80, 5), side=1, labels = TRUE)
  legend("topright", legend = c("absolute Best Set","within tolerance Loss "), pch = c(22,25), col=c("red","black"), pt.bg=c("red","black"), text.col=c("red","black"))
  
  ################
  ## Variable importance evaluation
  # Random Forest: from the R package
  selvarImp <- varImp(rfSelProfile, scale = TRUE)
  
  selfeatureswtol = data.frame(varImp = selvarImp[1:within5Pct,])
  selfeatureswtol$features = rownames(selvarImp)[1:within5Pct]
  
  output<-list(selfeatures=selfeatureswtol, selvarImp=selvarImp, within5Pct=within5Pct)
  return(output)  
}
