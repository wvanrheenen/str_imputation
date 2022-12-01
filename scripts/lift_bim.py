#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:43:01 2018

@author: wrheene2
"""

import sys
import pandas as pd
from pyliftover import LiftOver

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
    lo = LiftOver("hg19", "hg38")

    def lift(pos):
        pos_new = lo.convert_coordinate("chr" + str(pos[0]), int(pos[1]))
        if pos_new == None or len(pos_new) == 0:
            bp_new = "NA"
            chrom_new = "NA"
        else:
            chrom_new = pos_new[0][0].replace("chr", "")
            bp_new = pos_new[0][1]
        if chrom_new not in [str(x) for x in list(range(1,25))]:
            bp_new = "NA"
            chrom_new = "NA"        
        return tuple([chrom_new, bp_new])


    data = pd.read_csv(snakemake.input[0], sep="\s+", header=None)

    old_pos = list(zip(data[0], data[3]))
    new_pos = list(map(lift, old_pos))

    data[0], data[3] = zip(*new_pos)

    out_lifted = data.loc[ ~(data[0] == "NA")]
    out_notlifted = data.loc[ data[0] == "NA"][1]


    out_lifted.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
    out_notlifted.to_csv(snakemake.output[1], sep="\t", index=False, header=False)
