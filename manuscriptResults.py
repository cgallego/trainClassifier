# -*- coding: utf-8 -*-
"""
Classifier comparisons for 2 stage cascade vs. one-shot all benign all malignant

@author:  Cristina
"""
from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
from rpy2.robjects import globalenv
devtools = importr("devtools")
utils = importr('utils')
base = importr('base')

base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/compareRF_cascadeTest.R')
base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/display_cascadeTest.R')
base.source('Z:/Cristina/MassNonmass/codeProject/codeBase/trainClassifier/Rcode/ci_cascadeTest.R')

features = utils.read_table("4class_allFeatures.txt", header=True)
print(features.colnames)
print(base.nrow(features))

output = globalenv['compareRF_cascadeTest'](2, features, verbose=True)

results_cascade_vs_oneshot =  globalenv['display_cascadeTest'](output, Sen_cascade=False, Spec_cascade=False, printpvals=True)
print( results_cascade_vs_oneshot.rx2("pvaltable") )

# calculate CI results
globalenv['ci_cascadeTest'](output)

##### retrieve results
cascadeCollected = output.rx2("cascadeCollected")
pred1 = list( cascadeCollected[cascadeCollected.names.index('pred1')])