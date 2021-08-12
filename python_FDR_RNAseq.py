import codecs
import statsmodels.stats.multitest as smm
import csv
import numpy as np

pvalues = []
with codecs.open('/Users/SLancaster/Desktop/SCInulin_vs_Arabinoxylan_pvalues.csv', "r", encoding="utf-8-sig") as csvfile:
	spamreader = csv.reader(csvfile)
	for row in spamreader:
		pvalues.append(float(row[0]))

pvalues = np.array(pvalues)
pvalue_corr = smm.multipletests(pvalues, alpha = 0.05, method="fdr_tsbhy")[1] #correcting each individual pvalue using the benjamini hochberg correction
