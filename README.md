# LSM3241_CA1
R, Python scripts and supplementary information for LSM3241 Genomic Data Analysis CA1

## Pre-processing & Linear Models 
### linear-model-analysis.R 
This file allows analysis for an affymatrix microarray with linear models being built using the limma package. Gene ontology analysis is also included in this script. 

### preprocessing-of-data-for-meta-analysis.R 
This file reads in .CEL files for a particular GEO series and extracts gene expression data following RMA normalisation. 

## Gene-Ontology-Visualisation 
This file shows the output of the gene ontology based enrichment analysis for three nodes. 

## Meta-Analysis with Weka
### build-test-data
Build testing data from the original data set
#### convert-arff.py
Read the gse*.csv data set and output an ARFF file for weka, test_data.arff
#### gse50697.csv

### build-train-data
Build training data from the 5 data sets 
#### convert-arff.py
Read the gse*.csv dataset and output an ARFF file for weka, train_data.arff
#### gse12548.csv
#### gse14405.csv
#### gse14773.csv
#### gse18070.csv
#### gse43489.csv

## results
Results build using train_data.arff
Classifier used: Sequential Minimal Optimization (SMO) classifier, which is a a variant of Support Vector Machine (SVM) classifier, on WEKA 3.8.
The SMO classifier was built with a polynomial kernel (linear; e=1), a complexity parameter of 1 (c=1) and a logistic regression calibrator. The classification was done with 10-fold cross-validation.
#### SMO-classification-test-details.txt
#### SMO-classification-test-results.csv
#### SMO-classification-train-details.txt
#### SMO-classification-train-results.csv
#### features-weights-processing.py
Process feature weights from SMO and output the following:
Distribution of absolute feature weights
MLE-derived Exponential Distribution
Probability of observing a feature weight with the MLE-derived exponential distribution
