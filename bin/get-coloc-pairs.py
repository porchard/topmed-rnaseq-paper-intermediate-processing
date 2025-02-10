#!/usr/bin/env python
# coding: utf-8


import sys
import json
import pandas as pd
import pybedtools as bt


def get_pairs(first, second):
    first_bt = pd.DataFrame([[k] + v['region'].split('_') for k, v in first.items()], columns=['susie', 'chrom', 'start', 'end'])
    first_bt.start = first_bt.start.astype(int)
    first_bt.end = first_bt.end.astype(int)
    second_bt = pd.DataFrame([[k] + v['region'].split('_') for k, v in second.items()], columns=['susie', 'chrom', 'start', 'end'])
    second_bt.start = second_bt.start.astype(int)
    second_bt.end = second_bt.end.astype(int)
    first_bt = bt.BedTool.from_dataframe(first_bt[['chrom', 'start', 'end', 'susie']]).sort()
    second_bt = bt.BedTool.from_dataframe(second_bt[['chrom', 'start', 'end', 'susie']]).sort()
    overlaps = first_bt.intersect(second_bt, wa=True, wb=True)
    return overlaps.to_dataframe(names=['chrom', 'start', 'end', 'susie1', 'chrom2', 'start2', 'end2', 'susie2'])


first = json.load(open(sys.argv[1], 'r'))
second = json.load(open(sys.argv[2], 'r'))

tmp = get_pairs(first, second)

if len(tmp) > 0:
    tmp.loc[:,['susie1', 'susie2']].to_csv(sys.stdout, sep='\t', index=False, header=False)