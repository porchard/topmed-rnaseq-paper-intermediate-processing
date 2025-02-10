#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import sys

MATRICES = sys.argv[1:]


matrices = [pd.read_csv(f, sep='\t', index_col=0) for f in MATRICES]
column_counts = {}
for m in matrices:
    for c in m.columns:
        if c not in column_counts:
            column_counts[c] = 0
        column_counts[c] += 1
keep_columns = [c for c, count in column_counts.items() if count == len(matrices)]
combined = pd.concat([m[keep_columns] for m in matrices])
combined.to_csv(sys.stdout, sep='\t', index=True, header=True)