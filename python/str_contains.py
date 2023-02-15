import csv
from numpy import loadtxt
import pandas as pd
import numpy as np

fname = '/Users/sq566/Downloads/nsrr_studies/apples/data/edf/subjects.txt'
subjects = loadtxt(fname, dtype='str').tolist()

for sub in subjects:

    in_csv_file = f'/Users/sq566/Downloads/nsrr_studies/apples/data/csv/modified/{sub}.csv'
    df = pd.read_csv(in_csv_file, encoding= 'unicode_escape', header=None)
    
    d1 = df.iloc[:,0].str.contains(r'\d{1,2}-\d{1,2}-\d{2}', regex=True, na=True)
    
    d2 = df.iloc[:,0].str.contains(r'\d{1,2}/\d{1,2}/\d{2}', regex=True, na=True)
    
    if sum(d1) > 1 or sum(d2) > 1:
        print(sub)
        print(np.flatnonzero(d1))
        print(np.flatnonzero(d2))
        print("------------------------------------------------------------")
