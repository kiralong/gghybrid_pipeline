#Performing a Fisher's Exact Test on the Manacus hybrid zone cline enrichment data
#For the Candidate Gene Paper.
#5/19/2021

#The matrix needs to look like the following (with number examples from counting my r25 cline enrichment dataset for cline CENTERS):
#Significant cline loci in genome outlier windows[13]   sig. clines NOT in genome outlier windows[3,279]
#Non-significant clines in genome outlier windows[76]   Non-significant loci NOT in genome outlier windows[32,714]

#Make an object that contains the matrix of interest containing the above info
manacus_clines_r25_c <-  matrix(c(13,76,3279,32714), nrow = 2, ncol = 2)
manacus_clines_r25_c

#Now do the test!
FEmanacus_r25_c <- fisher.test(manacus_clines_r25_c)
FEmanacus_r25_c

#The results output
#Fisher's Exact Test for Count Data

#data:  manacus_clines_r25_c
#p-value = 0.09296
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#0.8681075 3.1002755
#sample estimates:
#odds ratio 
#1.706524 

#Now rinse and repeat to test the rest of the cline enrichment data

#Make an object that contains the matrix of interest for the r50 dataset on cline centers with pops 3&4
manacus_clines_r50_c <-  matrix(c(14,38,3095,15668), nrow = 2, ncol = 2)
manacus_clines_r50_c

#Now do the test!
FEmanacus_r50_c <- fisher.test(manacus_clines_r50_c)
FEmanacus_r50_c

#Fisher's Exact Test for Count Data
#
#data:  manacus_clines_r50_c
#p-value = 0.05865
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#0.9322109 3.5263601
#sample estimates:
#odds ratio 
#1.864967 

#Make an object that contains the matrix of interest for the r50 dataset on cline width with pops 3&4
manacus_clines_r50_w <-  matrix(c(19,33,4863,13900), nrow = 2, ncol = 2)
manacus_clines_r50_w

#Now do the test!
FEmanacus_r50_w <- fisher.test(manacus_clines_r50_w)
FEmanacus_r50_w

#Fisher's Exact Test for Count Data
#
#data:  manacus_clines_r50_w
#p-value = 0.08332
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#0.8832403 2.9834701
#sample estimates:
#odds ratio 
#1.645649 

