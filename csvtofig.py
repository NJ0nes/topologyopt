import csv
import matplotlib.pyplot as pyplot
import sys

table = []

with open(sys.argv[1], mode = "r") as f:
    in_data = csv.reader(f)
    table = []
    for line in in_data:
        table += [[float(s) for s in line]]

m = max([max(n) for n in table])

pyplot.imshow(table, cmap='gray_r', vmin=0.0, vmax=1.0)
pyplot.savefig(sys.argv[2])
