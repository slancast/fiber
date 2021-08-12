#! usr/bin/python3
"""
This class performed random forest training on the metabolomics dataset.
"""
#A simple find and replace for the fiber types works to iterate over the three of them
#i.e. Arabinoxylan, LCInulin, Mix

import csv
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import recall_score
from sklearn.model_selection import KFold
from collections import defaultdict
from collections import OrderedDict
import textwrap

import statistics

import scipy
import statistics
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


df = pd.read_csv('/Users/SLancaster/Desktop/Supervised_Clustering/New_Supclust_Datafiles/rna/rna_for_supervised_Arabinoxylan_binary.txt',sep='\t')
timepoint = np.array(df["binary"])
df = df.drop(['week','participant','binary'], axis = 1)
df = df.apply(lambda x: x.fillna(x.mean()),axis=0)

# Saving feature names for later use
headers_list = list(df.columns)
# Convert to numpy array
df = np.array(df)

kfold = KFold(5)
total_feature_importances = defaultdict(list)
median_feature_imporances = {} #It looks like Anna and I talked about median rather than mean. Maybe that is the way to go.
accuracies = []
recalls = []

for train_index, test_index in kfold.split(df):
    X_train, X_test = df[train_index], df[test_index]
    y_train, y_test = timepoint[train_index], timepoint[test_index]
    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(X_train, y_train)
    predictions = rf.predict(X_test)
    importances = list(rf.feature_importances_)
    feature_importances = [(feature, importance) for feature, importance in zip(headers_list, importances)]
    for i in feature_importances:
        total_feature_importances[str(i[0])].append(i[1])
    feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
    accuracy = accuracy_score(y_test, predictions)
    accuracies.append(accuracy)
    recall = recall_score(y_test, predictions)
    recalls.append(recall)

pd_feature_importances = pd.DataFrame.from_dict(total_feature_importances, orient='index')
importances_means = pd_feature_importances.mean(axis=1)
importances_means.to_csv("/Users/SLancaster/Desktop/rna_Arabinoxylan_Importances.csv")

accuracies_recalls = pd.DataFrame(accuracies, index=["cross_validation_1","2","3","4","5"], columns = ["Accuracies"])
accuracies_recalls['Recalls'] = recalls
accuracies_recalls.to_csv("/Users/SLancaster/Desktop/rna_Arabinoxylan_Accuracies_Recalls.csv")

for i in total_feature_importances:
    median_feature_imporances[i] = statistics.median(total_feature_importances[i])

ordered_importances = OrderedDict(sorted(median_feature_imporances.items(), key=lambda t: t[1], reverse=True))
top_20_features = []
counter = 0
for k, v in ordered_importances.items():
    if counter <= 20:
        top_20_features.append(k)
    counter += 1


numpy_array = []
for i in range(0, 21):
    numpy_array.append(total_feature_importances[top_20_features[i]])

numpy_array = np.array(numpy_array)

nums = ["1", "2", "3", "4", "5"]

fig, ax = plt.subplots()
im = ax.imshow(numpy_array)

# We want to show all ticks...
ax.set_xticks(np.arange(len(nums)))
ax.set_yticks(np.arange(len(top_20_features)))
# ... and label them with the respective list entries
ax.set_xticklabels(nums)
ax.set_yticklabels(top_20_features, {'fontsize': 6})

# Create colorbar
cbar = ax.figure.colorbar(im, ax=ax)
#cbar.ax.set_ylabel(rotation=-90, va="bottom")

# Rotate the tick labels and set their alignment.
plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
         rotation_mode="anchor")

# Loop over data dimensions and create text annotations.


#plt.figure(figsize=(10,5))
ax.set_title("Feature's Impact on rna_Arabinoxylan")
fig.tight_layout()
plt.show()
fig.savefig("/Users/SLancaster/Desktop/rna_Arabinoxylan.png")

