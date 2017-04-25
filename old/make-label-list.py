#! /usr/bin/env python
from __future__ import print_function
import sys, screed

label_list = []
for record in screed.open(sys.argv[1]):
    label_list.append(record.name)

for n, l in enumerate(label_list):
    print("{} {}\n".format(n + 1, l))
