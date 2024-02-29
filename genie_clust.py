import genieclust
import numpy as np
import matplotlib.pyplot as plt
import csv

csv_filename = "outcome.csv"
csv_filename2 = "ground_truth_index16.csv"


a = []

with open(csv_filename, 'r', newline='') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        if len(row) >= 3:
            a.append(int(float(row[2].strip())))
            # a.append(int(row[2].strip()))


b = []
with open(csv_filename2, 'r', newline='') as csv_file:
    csv_reader2 = csv.reader(csv_file)
    for row in csv_reader2:
        b.append(row[0])


# g = genieclust.Genie(n_clusters=5)
# labels = g.fit_predict(points)
# genieclust.plots.plot_scatter(points, labels=labels, alpha=0.5)
# genieclust.compare_partitions.normalized_accuracy

# a=[1,1,2,1,2,3,1,3,2,1]
# b=[1,1,2,1,2,3,1,3,2,1]
# print(a)
# print(b)

print(genieclust.compare_partitions.normalized_clustering_accuracy(a, b))

# plt.show()
