import pandas as pd
from numpy import loadtxt

fname = '/Users/sq566/Desktop/mesa/subjects.txt'
subjects = loadtxt(fname, dtype='str').tolist()

data = {'class': [], 'instance': [], 'channel': [], 'start': [], 'stop': [], 'meta': []}

for ID in subjects:
    
    file_path =  f'/Users/sq566/Desktop/mesa/xls/{ID}.xls'
    print("Processing...", file_path)
    
    df = pd.read_excel(file_path, engine='xlrd', usecols=['Event', 'Start Time', 'End Time'])
    df = df.iloc[1:].reset_index(drop=True)
    
    df['Start Time'] = df['Start Time'].astype(str)
    df['End Time'] = df['End Time'].astype(str)
    df['Start Time'] = df['Start Time'].apply(lambda x: f"{x[8:10]}-{x[5:7]}-{x[0:4]}-{x[11:]}")
    df['End Time'] = df['End Time'].apply(lambda x: f"{x[8:10]}-{x[5:7]}-{x[0:4]}-{x[11:]}")
    
    
    df_out = pd.DataFrame(data)
    df_out['class'] = df['Event']
    df_out['instance'] = "."
    df_out['channel'] = "."
    df_out['start'] = df['Start Time']
    df_out['stop'] = df['End Time']
    df_out['meta'] = "."
    
    df_out.to_csv(f'/Users/sq566/Desktop/mesa/annots/{ID}.annot', index=None, sep='\t', header=None)
    del df_out, df
