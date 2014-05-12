display_cascadeTest <- function(outputresamples, Sen_cascade=FALSE, Spec_cascade=FALSE, printpvals=TRUE, multi=FALSE) {
  ####### AUC
  RF1AUCframe=data.frame(test = 1:length(outputresamples$RFmodel1_AUC))
  RF1AUCframe$test = (outputresamples$RFmodel1_AUC)/100
  RF1AUCframe$MorN = "stage1"
  RF1AUCframe$metric = "AUC"
  RF1AUCframe$model = "RF"
  allmetricsframe = RF1AUCframe
  print("stage1 AUC")
  print( median(outputresamples$RFmodel1_AUC/100) ) 
  print( quantile(outputresamples$RFmodel1_AUC/100) ) 
  
  RF2AUCframe=data.frame(test = 1:length(outputresamples$RFmodel2_AUC))
  RF2AUCframe$test = (outputresamples$RFmodel2_AUC)/100
  RF2AUCframe$MorN = "stage2_mass"
  RF2AUCframe$metric = "AUC"
  RF2AUCframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF2AUCframe),] <-  RF2AUCframe
  print("stage2_mass AUC")
  print( median(outputresamples$RFmodel2_AUC/100) ) 
  print( quantile(outputresamples$RFmodel2_AUC/100) ) 
  
  RF3AUCframe=data.frame(test = 1:length(outputresamples$RFmodel3_AUC))
  RF3AUCframe$test = (outputresamples$RFmodel3_AUC)/100
  RF3AUCframe$MorN = "stage2_nonmass"
  RF3AUCframe$metric = "AUC"
  RF3AUCframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF3AUCframe),] <-  RF3AUCframe
  print("stage2_nonmass AUC")
  print( median(outputresamples$RFmodel3_AUC/100) ) 
  print( quantile(outputresamples$RFmodel3_AUC/100) ) 
  
  ####### Sen
  RF1Senframe=data.frame(test = 1:length(outputresamples$RFmodel1_Sen))
  RF1Senframe$test = outputresamples$RFmodel1_Sen
  RF1Senframe$MorN = "stage1"
  RF1Senframe$metric = "Sen"
  RF1Senframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF1Senframe),] <-  RF1Senframe
  print("stage1 Sen")
  print( median(outputresamples$RFmodel1_Sen) )
  print( quantile(outputresamples$RFmodel1_Sen) )
  
  RF2Senframe=data.frame(test = 1:length(outputresamples$RFmodel2_Sen))
  RF2Senframe$test = outputresamples$RFmodel2_Sen
  RF2Senframe$MorN = "stage2_mass"
  RF2Senframe$metric = "Sen"
  RF2Senframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF2Senframe),] <-  RF2Senframe
  print("stage2_mass Sen")
  print( median(outputresamples$RFmodel2_Sen) )
  print( quantile(outputresamples$RFmodel2_Sen) )
  
  RF3Senframe=data.frame(test = 1:length(outputresamples$RFmodel3_Sen))
  RF3Senframe$test = outputresamples$RFmodel3_Sen
  RF3Senframe$MorN = "stage2_nonmass"
  RF3Senframe$metric = "Sen"
  RF3Senframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF3Senframe),] <-  RF3Senframe
  print("stage2_nonmass Sen")
  print( median(outputresamples$RFmodel3_Sen) )
  print( quantile(outputresamples$RFmodel3_Sen) )
  
  ####### Spec
  RF1Specframe=data.frame(test = 1:length(outputresamples$RFmodel1_Spec))
  RF1Specframe$test = outputresamples$RFmodel1_Spec
  RF1Specframe$MorN = "stage1"
  RF1Specframe$metric = "Spec"
  RF1Specframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF1Specframe),] <-  RF1Specframe
  print("stage1 Spec")
  print( median(outputresamples$RFmodel1_Spec) )
  print( quantile(outputresamples$RFmodel1_Spec) )
  
  RF2Specframe=data.frame(test = 1:length(outputresamples$RFmodel2_Spec))
  RF2Specframe$test = outputresamples$RFmodel2_Spec
  RF2Specframe$MorN = "stage2_mass"
  RF2Specframe$metric = "Spec"
  RF2Specframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF2Specframe),] <-  RF2Specframe
  print("stage2_mass Spec")
  print( median(outputresamples$RFmodel2_Spec) )
  print( quantile(outputresamples$RFmodel2_Spec) )
  
  RF3Specframe=data.frame(test = 1:length(outputresamples$RFmodel3_Spec))
  RF3Specframe$test = outputresamples$RFmodel3_Spec
  RF3Specframe$MorN = "stage2_nonmass"
  RF3Specframe$metric = "Spec"
  RF3Specframe$model = "RF"
  allmetricsframe[nrow(allmetricsframe)+1:nrow(RF3Specframe),] <-  RF3Specframe
  print("stage2_nonmass Spec")
  print( median(outputresamples$RFmodel3_Spec) )
  print( quantile(outputresamples$RFmodel3_Spec) )
  
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  
  print(
    ggplot(data = allmetricsframe, aes(x = metric, y = test)) +
      geom_boxplot(aes(fill = interaction(model, MorN),
                       group = interaction(factor(metric), model, MorN)))+
      theme_bw(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(name = "2-stage individual classifier performance",
                        values = c("#99CC00","#993366","#008080")) 
  )
  
  
  # Now plot comparisons
  ###########
  ####### AUC
  RFcascadeAUCframe=data.frame(test = 1:length(outputresamples$cascadeRFmodels_AUC))
  RFcascadeAUCframe$test = (outputresamples$cascadeRFmodels_AUC)/100
  RFcascadeAUCframe$MorN = "cascade"
  RFcascadeAUCframe$metric = "AUC"
  RFcascadeAUCframe$model = "RF"
  allmetricsframe_cascade = RFcascadeAUCframe
  print("cascade AUC")
  print( median(outputresamples$cascadeRFmodels_AUC/100) )
  print( quantile(outputresamples$cascadeRFmodels_AUC/100) )
  
  RFoneshotAUCframe=data.frame(test = 1:length(outputresamples$oneshotRFmodel4_AUC))
  RFoneshotAUCframe$test = (outputresamples$oneshotRFmodel4_AUC)/100
  RFoneshotAUCframe$MorN = "oneshot"
  RFoneshotAUCframe$metric = "AUC"
  RFoneshotAUCframe$model = "RF"
  allmetricsframe_cascade[nrow(allmetricsframe_cascade)+1:nrow(RFoneshotAUCframe),] <-  RFoneshotAUCframe
  print("oneshot AUC")
  print( median(outputresamples$oneshotRFmodel4_AUC/100) )
  print( quantile(outputresamples$oneshotRFmodel4_AUC/100) )
  
  ####### Sen
  if(Sen_cascade){
    RFcascadeSenframe=data.frame(test = 1:length(outputresamples$cascadeRFmodels_Sen))
    RFcascadeSenframe$test = outputresamples$cascadeRFmodels_Sen
    RFcascadeSenframe$MorN = "cascade"
    RFcascadeSenframe$metric = "Sen"
    RFcascadeSenframe$model = "RF"
    allmetricsframe_cascade[nrow(allmetricsframe_cascade)+1:nrow(RFcascadeSenframe),] <-  RFcascadeSenframe
    print("cascade Sen")
    print( median(outputresamples$cascadeRFmodels_Sen) )
    print( quantile(outputresamples$cascadeRFmodels_Sen) )
    
    RFoneshotSenframe=data.frame(test = 1:length(outputresamples$oneshotRFmodel4_Sen))
    RFoneshotSenframe$test = outputresamples$oneshotRFmodel4_Sen
    RFoneshotSenframe$MorN = "oneshot"
    RFoneshotSenframe$metric = "Sen"
    RFoneshotSenframe$model = "RF"
    allmetricsframe_cascade[nrow(allmetricsframe_cascade)+1:nrow(RFoneshotSenframe),] <-  RFoneshotSenframe
    print("oneshot Sen")
    print( median(outputresamples$oneshotRFmodel4_Sen) )
    print( quantile(outputresamples$oneshotRFmodel4_Sen) )
  }
  ####### Spec
  if(Spec_cascade){
    RFcascadeSpecframe=data.frame(test = 1:length(outputresamples$cascadeRFmodels_Spec))
    RFcascadeSpecframe$test = outputresamples$cascadeRFmodels_Spec
    RFcascadeSpecframe$MorN = "cascade"
    RFcascadeSpecframe$metric = "Spec"
    RFcascadeSpecframe$model = "RF"
    allmetricsframe_cascade[nrow(allmetricsframe_cascade)+1:nrow(RFcascadeSpecframe),] <-  RFcascadeSpecframe
    print("cascade Spec")
    print( median(outputresamples$cascadeRFmodels_Spec) )
    print( quantile(outputresamples$cascadeRFmodels_Spec) )
    
    RFoneshotSpecframe=data.frame(test = 1:length(outputresamples$oneshotRFmodel4_Spec))
    RFoneshotSpecframe$test = outputresamples$oneshotRFmodel4_Spec
    RFoneshotSpecframe$MorN = "oneshot"
    RFoneshotSpecframe$metric = "Spec"
    RFoneshotSpecframe$model = "RF"
    allmetricsframe_cascade[nrow(allmetricsframe_cascade)+1:nrow(RFoneshotSpecframe),] <-  RFoneshotSpecframe
    print("oneshot Spec")
    print( median(outputresamples$oneshotRFmodel4_Spec) )
    print( quantile(outputresamples$oneshotRFmodel4_Spec) )
  }
  
  
  print(
    ggplot(data = allmetricsframe_cascade, aes(x = metric, y = test)) +
      geom_boxplot(aes(fill = interaction(model, MorN),
                       group = interaction(factor(metric), model, MorN)))+
      theme_bw(base_size = 16) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_manual(name = "Cascade vs. Oneshot Results",
                        values = c("#99CC00","#993366","#008080")) )
  ##############
  ### Calculate resample statistics for cascade vs. one shot
  if(printpvals){
    RFAUC_ttest = t.test((RFcascadeAUCframe$test)*100, (RFoneshotAUCframe$test)*100, alternative = "greater", paired=TRUE, var.equal = FALSE)
    print(RFAUC_ttest)
    pvalframe=data.frame(pvalue=RFAUC_ttest$p.value)
    pvalframe$test = "AUC"
    
    if(Sen_cascade){
      RFSen_ttest = t.test((RFcascadeSenframe$test)*100, (RFoneshotSenframe$test)*100, alternative = "greater", paired=TRUE, var.equal = FALSE)
      print(RFSen_ttest)
      pvalframe[2,1]=RFSen_ttest$p.value
      pvalframe[2,2] = "Sen"
    }
    if(Spec_cascade){
      RFSpec_ttest = t.test((RF1Specframe$test)*100, (RF2Specframe$test)*100, alternative = "greater", paired=FALSE, var.equal = FALSE)
      print(RFSpec_ttest)
      pvalframe[3,1]=RFSpec_ttest$p.value
      pvalframe[3,2] = "Spec"
    }
    
    pval <- list(pvalue = paste0("pval: ", format.pval(pvalframe[,1], eps=0.001)))
    pvaltable=data.frame(pvalue=pval)
  }
  
  if(printpvals && Sen_cascade && Spec_cascade){
    print(
      ggplot(data = allmetricsframe_cascade, aes(x = metric, y = test)) +
        geom_boxplot(aes(fill = interaction(model, MorN),
                         group = interaction(factor(metric), model, MorN)))+
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(name = "Cascade vs. Oneshot Results",
                          values = c("#99CC00","#993366","#008080")) +
        annotate("text", label = pvaltable[1,1], x = 1, y = 1, size = 6) +
        annotate("text", label = pvaltable[2,1], x = 2, y = 1, size = 6) +
        annotate("text", label = pvaltable[3,1], x = 3, y = 1, size = 6)
    )
  }
  if(printpvals && Sen_cascade==FALSE && Spec_cascade==FALSE ){
    print(
      ggplot(data = allmetricsframe_cascade, aes(x = metric, y = test)) +
        geom_boxplot(aes(fill = interaction(model, MorN),
                         group = interaction(factor(metric), model, MorN)))+
        theme_bw(base_size = 16) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_manual(name = "Cascade vs. Oneshot Results",
                          values = c("#99CC00","#993366","#008080")) +
        annotate("text", label = pvaltable[1,1], x = 1, y = 1, size = 6)
    )
  }
  
  output<-list(allmetricsframe=allmetricsframe, allmetricsframe_cascade=allmetricsframe_cascade, pvaltable=pvaltable)
  
  return(output)  
}
