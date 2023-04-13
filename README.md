# AML-KNN/NSC
Subtype classification of Acute Myeloid Leukemia samples using k nearest neighbors with nearest shrunken centroids for data reduction

This script will take an AML subtype as input and return a dataframe df.results.subtype with 
the results from a 5x5-fold cross validated NSC job, to create a classifier for the specified 
subtype. It will also return a list of cpg sites with their frequency of occurrence in the 25 
runs. Specify the AML subtype under "Preparing to run the script".

The results from the script will be stored in two variables: df.results.subtype and cpg.freq.subtype
df.results will have 25 rows (one for each run) and 6 columns:
"threshold" is the threshold used for NSC shrinking of centroids
"N_CpGs" is the number of CpG sites used to classify the subtype using the chosen threshold
"error_rate_conf" is the number of errors seen in the confusion table when predicting class
"error_predict" is the number of errors that arise when testing the classifier on the test set
"CpG_IDs" is a list of indices for the chosen CpG sites
"test_set" is a list of indices for the samples that are included in the test set for the run

cpg.freq.subtype is an ordered dataframe with the indices of all chosen CpG sites for all 25 runs
and the frequency with which they were chosen. 25 would been chosen in all runs, 8 would mean
chosen in 8 runs out of 25.
