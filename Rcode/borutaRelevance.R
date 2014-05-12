findRelevant <- function(massallfeatures, nonmassallfeatures) {
    library("Boruta")
    require(data.table)
    require(ggplot2)
        
    #Boruta mass:
    set.seed(111)
    massBoruta <-Boruta(BenignNMaligNAnt ~., data=massallfeatures, doTrace=2, ntree=500)
    massBoruta
    
    plot(massBoruta)
    #Color codes:
    #  c("green", "yellow", "red", "blue"),
    #  Confirmed, Tentative, Rejected and shadow.
    # Blue boxplots correspond to minimal, average and maximum Z score of a shadow attribute. Red
    
    # Short
    shortmassBoruta <-Boruta(BenignNMaligNAnt ~., data=massallfeatures, doTrace=2, maxRuns=12)
    plot.new()
    plot(shortmassBoruta)
    
    #Boruta non-mass:
    set.seed(111)
    nonmassBoruta <-Boruta(BenignNMaligNAnt ~., data=nonmassallfeatures, doTrace=2, ntree=500)
    nonmassBoruta
    
    plot.new()
    plot(nonmassBoruta)
    
    shortmassBoruta <-Boruta(BenignNMaligNAnt ~., data=massallfeatures, doTrace=2, maxRuns=12)
    plot.new()
    plot(shortmassBoruta)
    
    ## Now format data and collect relevant features Z-scores:
    
    library("latticeExtra")
    rankingsmass <- massBoruta$ImpHistory
    rankingsnonmass <- nonmassBoruta$ImpHistory
    
    confirmedmass_features<-massBoruta$finalDecision[massBoruta$finalDecision == "Confirmed"]
    # Confirmed mass features
    print(confirmedmass_features)
    
    confirmednonmass_features<-nonmassBoruta$finalDecision[nonmassBoruta$finalDecision == "Confirmed"]
    # Confirmed nonmass features
    print(confirmednonmass_features)
    
    ####### proces Masses (add fist confirmed feature)
    cfeature = as.data.frame(confirmedmass_features[1])  
    massframe=data.frame(zscores=rankingsmass[is.finite(rankingsmass[,rownames(cfeature)]),rownames(cfeature)])
    massframe$MorN = "mass"
    massframe$feature = rownames(cfeature)
    masszscore_selected <- massframe
    
    nonmassframe=data.frame(zscores=rankingsnonmass[is.finite(rankingsnonmass[,rownames(cfeature)]),rownames(cfeature)])
    nonmassframe$MorN = "nonmass"
    nonmassframe$feature = rownames(cfeature)
    masszscore_selected[nrow(masszscore_selected)+1:nrow(nonmassframe),] <-  nonmassframe
    
    masszscore_ttest = numeric(length(confirmedmass_features))
    masszscore_ttest[1] = t.test(as.data.frame(massframe)$zscores, as.data.frame(nonmassframe)$zscores, alternative = "greater", paired=FALSE, var.equal = FALSE)$p.value
    
    pvallabels = character(length(confirmedmass_features))
    pvallabels[1] <-  rownames(cfeature)
    
    # proces remaining confirmed feature Masses
    for (i in 2:length(confirmedmass_features)) {
      cfeature = as.data.frame(confirmedmass_features[i])  
      massframe=data.frame(zscores=rankingsmass[,rownames(cfeature)])
      massframe$MorN = "mass"
      massframe$feature = rownames(cfeature)
      masszscore_selected[nrow(masszscore_selected)+1:nrow(massframe),] <-  massframe
      
      nonmassframe=data.frame(zscores=rankingsnonmass[,rownames(cfeature)])
      nonmassframe$MorN = "nonmass"
      nonmassframe$feature = rownames(cfeature)
      masszscore_selected[nrow(masszscore_selected)+1:nrow(nonmassframe),] <-  nonmassframe
      
      # p-values test
      masszscore_ttest[i] = t.test(massframe$zscores[is.finite(massframe$zscores)], nonmassframe$zscores[is.finite(nonmassframe$zscores)], alternative = "greater", paired=FALSE, var.equal = FALSE)$p.value
      pvallabels[i] <-  rownames(cfeature)
    }
    
    #Now plot the compasion with respect to non-mass features for relevant mass features:
    # Done processing
    print(confirmedmass_features)
    
    # format p values and order in dataframe
    dt <- data.table(masszscore_ttest)
    pval <- dt[, list(pvalue = paste0("pval: ", format.pval(masszscore_ttest,  eps=0.001)))]
    
    pvalframe=data.frame(pvalue=pval)
    pvalframe$feature = pvallabels
    order=order(pvalframe$pvalue, pvalframe$feature)
    order[5]=4
    order[6]=5
    pvalframe[order,]
    
    # generate boxplot comparison of relevant mass features vs. the same non-mass feature
    plot.new()
    ggplot() + 
      geom_boxplot(data = masszscore_selected, mapping=aes(x=MorN, y=zscores, fill = factor(MorN))) + 
      geom_text(data = pvalframe, aes(label=pvalue, x="mass", y=9)) +
      facet_grid(~ feature) + theme_bw(base_size = 16) + 
      labs(title = "Comparison of Z-scores for Mass confirmed features", y="Z-scores")+
      theme(  axis.text.x=element_text(angle=0, face="bold", size=12),
              legend.position = "bottom",
              strip.background = element_rect(fill=NA), 
              panel.grid.major = element_line(colour=NA), 
              panel.grid.minor = element_line(colour=NA))
    
 
    output<-list(masszscore_selected=masszscore_selected, pvalframe=pvalframe)
    
    return(output)  
}