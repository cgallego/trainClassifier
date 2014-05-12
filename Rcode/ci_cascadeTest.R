ci_cascadeTest <-  function(outputresamples){
  ############### for AUC
  # append results for n=1500 bootstrap resamples
  MvsnonM_AUCmedianCIboot = apply(matrix(sample(outputresamples$RFmodel1_AUC, rep=TRUE, 10^4*length(outputresamples$RFmodel1_AUC)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel1_AUC))
  print(quantile(MvsnonM_AUCmedianCIboot, c(.025, 0.975)))
  
  onlymass_AUCmedianCIboot = apply(matrix(sample(outputresamples$RFmodel2_AUC, rep=TRUE, 10^4*length(outputresamples$RFmodel1_AUC)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel1_AUC))
  print(quantile(onlymass_AUCmedianCIboot, c(.025, 0.975)))
  
  onlynonmass_AUCmedianCIboot = apply(matrix(sample(outputresamples$RFmodel3_AUC, rep=TRUE, 10^4*length(outputresamples$RFmodel3_AUC)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel3_AUC))
  print(quantile(onlynonmass_AUCmedianCIboot, c(.025, 0.975)))
  
  
  ############### for Sensivity
  MvsnonM_SenmedianCIboot = apply(matrix(sample(outputresamples$RFmodel1_Sen, rep=TRUE, 10^4*length(outputresamples$RFmodel1_Sen)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel1_Sen))
  print(quantile(MvsnonM_SenmedianCIboot, c(.025, 0.975)))
  
  onlymass_SenmedianCIboot = apply(matrix(sample(outputresamples$RFmodel2_Sen, rep=TRUE, 10^4*length(outputresamples$RFmodel2_Sen)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel2_Sen))
  print(quantile(onlymass_SenmedianCIboot, c(.025, 0.975)))
  
  onlynonmass_SenmedianCIboot = apply(matrix(sample(outputresamples$RFmodel3_Sen, rep=TRUE, 10^4*length(outputresamples$RFmodel3_Sen)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel3_Sen))
  print(quantile(onlynonmass_SenmedianCIboot, c(.025, 0.975)))
  
  
  ############### for Specificity
  MvsnonM_SpecmedianCIboot = apply(matrix(sample(outputresamples$RFmodel1_Spec, rep=TRUE, 10^4*length(outputresamples$RFmodel1_Spec)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel1_Spec))
  print(quantile(MvsnonM_SpecmedianCIboot, c(.025, 0.975)))
  
  onlymass_SpecmedianCIboot = apply(matrix(sample(outputresamples$RFmodel2_Spec, rep=TRUE, 10^4*length(outputresamples$RFmodel2_Spec)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel2_Spec))
  print(quantile(onlymass_SpecmedianCIboot, c(.025, 0.975)))
  
  onlynonmass_SpecmedianCIboot = apply(matrix(sample(outputresamples$RFmodel3_Spec, rep=TRUE, 10^4*length(outputresamples$RFmodel3_Spec)), nrow=10^4), 1, median)
  print(summary(outputresamples$RFmodel3_Spec))
  print(quantile(onlynonmass_SpecmedianCIboot, c(.025, 0.975)))
  
  
  ######### compare AUC and deltaAUC cascade vs. oneshot ############
  cascade_AUCmedianCIboot = apply(matrix(sample(outputresamples$cascadeRFmodels_AUC, rep=TRUE, 10^4*length(outputresamples$cascadeRFmodels_AUC)), nrow=10^4), 1, median)
  print(summary(outputresamples$cascadeRFmodels_AUC))
  print(quantile(cascade_AUCmedianCIboot, c(.025, 0.975)))
  
  oneshot_AUCmedianCIboot = apply(matrix(sample(outputresamples$oneshotRFmodel4_AUC, rep=TRUE, 10^4*length(outputresamples$oneshotRFmodel4_AUC)), nrow=10^4), 1, median)
  print(summary(outputresamples$oneshotRFmodel4_AUC))
  print(quantile(oneshot_AUCmedianCIboot, c(.025, 0.975)))
  
  # now confidence intervals
  deltaAUC_medianCIboot = apply(matrix(sample(outputresamples$deltaAUC_RFcascadevsRFoneshot, rep=TRUE, 10^4*length(outputresamples$deltaAUC_RFcascadevsRFoneshot)), nrow=10^4), 1, median)
  print(summary(outputresamples$deltaAUC_RFcascadevsRFoneshot))
  print(quantile(deltaAUC_medianCIboot, c(.025, 0.975)))
  
  return()
}