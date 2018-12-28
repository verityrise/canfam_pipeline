#!/usr/bin/python
'''
This program is for integrating genome using glob function
Author: HZL
'''

import pandas as pd
import numpy as np
import matplotlib as mat
import os
import gzip
import glob
from integrate_genome import path, prefix

def qcplot(result):
    files=glob.glob(path+'*_CanFam3.1_*.fa.gz')
    outfile=gzip.open(path+prefix+".fa.gz",'wb')
    for f in files:
        gfile=gzip.open(f)
        outfile.write(gfile.read())
        gfile.close()


if __name__=="__main__":
    qcplot()