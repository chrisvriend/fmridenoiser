#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 12:52:28 2021

@author: C. Vriend - Amsterdam UMC 2022 
"""

import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser("determine number of motion outliers")
parser.add_argument("-dir", action="store", dest="dir", help="path to fmriprep directory", required=True)
parser.add_argument("-subjid", action="store", dest="subjid", help="subject ID", required=True)
parser.add_argument("-session", action="store", dest="session", help="session ID; if any", required=False)


if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()
fmriprepdir=args.dir
subj=args.subjid


if args.session is not None:
    session=args.session
    confoundfile=os.path.join(fmriprepdir,subj,session,'func',(subj + '_' + session + '_task-rest_desc-confounds_timeseries.tsv'))
    if not os.path.exists(confoundfile):
        print('cannot find confound tsv file')
        print(confoundfile)
        sys.exit(1)
    cols=pd.read_csv(confoundfile,sep='\t',nrows=1,header=None).transpose()
    outliers=cols[cols[0].str.match('motion_outlier')]
    outliers.to_csv(os.path.join(fmriprepdir,subj,session,'func',(subj + '_' + session + '_motion_outliers.txt')),index=False,header=False)
    print('extracted motion outlier columns to txt file')
    
    
else:
    confoundfile=os.path.join(fmriprepdir,subj,'func',subj + '_task-rest_desc-confounds_timeseries.tsv')
    if not os.path.exists(confoundfile):
        print('cannot find confound tsv file')
        print(confoundfile)
        sys.exit(1)
    cols=pd.read_csv(confoundfile,sep='\t',nrows=1,header=None).transpose()
    outliers=cols[cols[0].str.match('motion_outlier')]
    outliers.to_csv(os.path.join(fmriprepdir,subj,'func',(subj + '_motion_outliers.txt')),index=False,header=False)
    print('extracted motion outlier columns to txt file')

    


