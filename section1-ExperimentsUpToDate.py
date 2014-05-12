# -*- coding: utf-8 -*-
"""
Created on Tue May 06 15:30:56 2014

@ author (C) Cristina Gallego, University of Toronto
"""

import os, os.path
import sys
import string
from rpy2.robjects.packages import importr
import rpy2.robjects as R
from rpy2.robjects import globalenv
import rpy2.robjects.lib.ggplot2 as ggplot2

# convertion packages
import pandas.rpy.common as com
from rpy2.robjects.numpy2ri import numpy2ri


devtools = importr("devtools")
utils = importr('utils')
base = importr('base')
latt = importr('latticeExtra')
stats = importr('stats')


def runBoruta():
    base.load("Rcode/zscores.RData")
    base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/borutaRelevance.R')
    outputBoruta = globalenv['findRelevant'](globalenv['massallfeatures'], globalenv['nonmassallfeatures'])

    # generate boxplot comparison of relevant mass features vs. the same non-mass feature
    plotgp = ggplot2.ggplot(outputBoruta.rx2("masszscore_selected")) + \
          ggplot2.aes_string(x='MorN', y='zscores', fill = 'factor(MorN)') + \
          ggplot2.geom_boxplot() + \
          ggplot2.opts(title = "Comparison of Z-scores for Mass confirmed features", y="Z-scores") 
    plotgp.plot()
    
    return
    
def compareplotVarImp(selvarImp):
    base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/plotRFvarImp_2class.R')
    globalenv['plotRFvarImp_2class'](selvarImp)

    return    
    
def runFeatureSelection(flagOnly):
    base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/subset_feature_selection.R')
    
    # 1) onlyMass  
    #================
    if( flagOnly=="onlyMass" ):
        onlyMassfile="2classmass_allFeatures.txt"
        output_onlyMass = globalenv['subset_feature_selection'](onlyMassfile, 10)
        print( base.summary(output_onlyMass.rx2("selfeatures")))
                               
        selvarImpMass = output_onlyMass.rx2("selfeatures")
        output = selvarImpMass
    
    # 2) onlyNonmass
    #=================
    if( flagOnly=="onlyNonMass" ):
        onlyNonmassfile="2classnonmass_allFeatures.txt"
        output_onlyNonmass = globalenv['subset_feature_selection'](onlyNonmassfile, 1)
        print( base.summary( output_onlyNonmass.rx2("selfeatures")) )
        
        selvarImpNonMass = output_onlyNonmass.rx2("selfeatures")
        output = selvarImpNonMass
    
    # 3) Multi-class (4class):
    #=================
    if( flagOnly=="multiclass" ):
        multiclassfile="4class_allFeatures.txt"
        output_multiclass = globalenv['subset_feature_selection'](multiclassfile, 10)
        print( base.summary( output_multiclass.rx2("selfeatures")) )
        
        selvarImpmulticlass = output_multiclass.rx2("selfeatures")
        output = selvarImpmulticlass
        
    # 4) allBenignvsMalignant
    #=================
    if( flagOnly=="allBenignvsMalignant" ):
        allBenignvsMalignantfile="2class_allFeatures.txt"
        output_allBenignvsMalignant = globalenv['subset_feature_selection'](allBenignvsMalignantfile, 10)
        print( base.summary( output_allBenignvsMalignant.rx2("selfeatures")) )
        
        selvarImpallBenignvsMalignant = output_allBenignvsMalignant.rx2("selfeatures")
        output = selvarImpallBenignvsMalignant
        
    # 5) MassvsNonmass
    #=================
    if( flagOnly=="MassvsNonmass" ):
        MassvsNonmassfile="2class_RFvarImp_MorN.txt"
        output_MassvsNonmass  = globalenv['subset_feature_selection'](MassvsNonmassfile, 10)
        print( base.summary( output_MassvsNonmass.rx2("selfeatures")) )
        
        selvarImpMassvsNonmass = output_MassvsNonmass.rx2("selfeatures")
        output = selvarImpMassvsNonmass
    
    return output
    
    
    
if __name__ == '__main__':
    # Get Root folder ( the directory of the script being run)
    path_rootFolder = os.path.dirname(os.path.abspath(__file__))
    print path_rootFolder
    
    output_onlyMass = runFeatureSelection("onlyMass")
    output_onlyNonMass = runFeatureSelection("onlyNonMass")
    #selvarImpmulticlass = runFeatureSelection("multiclass")
    #selvarImpallBenignvsMalignant = runFeatureSelection("allBenignvsMalignant")
    #selvarImpMassvsNonmass = runFeatureSelection("MassvsNonmass")

    # Process for visualization
    pd_selvarImpMass = com.convert_robj(output_onlyMass)
    pd_selvarImpMass['SelectedFeatureGroup'] = "mass"
    selvarImpMass = com.convert_to_r_dataframe(pd_selvarImpMass)
    
    pd_selvarImpNonMass = com.convert_robj(output_onlyNonMass)
    pd_selvarImpNonMass['SelectedFeatureGroup'] = "nonmass"
    selvarImpNonMass = com.convert_to_r_dataframe(pd_selvarImpNonMass)
    allselvarImp = base.rbind(selvarImpMass, selvarImpNonMass)
    
    # DO a varImp plot comparison    
    compareplotVarImp(allselvarImp)
    
