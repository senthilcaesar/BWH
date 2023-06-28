import pandas as pd
import re
from numpy import loadtxt


fname = '/Users/sq566/Desktop/wsc/logs/logs.txt'
subjects = loadtxt(fname, dtype='str').tolist()

searchfor1 = ['LIGHTS OUT', 'Light Out']
searchfor2 = ['LIGHTS ON', 'Lights On']
words1 = [rf'\b{string}\b' for string in searchfor1]
words2 = [rf'\b{string}\b' for string in searchfor2]

df_out = pd.DataFrame(index=range(0,len(subjects)), columns=['ID','LIGHTS_OFF','LIGHTS_ON'])
    
for i, ID in enumerate(subjects):

    df_out.loc[i]['ID'] = ID
    filename = f'/Users/sq566/Desktop/wsc/logs/{ID}.log.txt'
    df = pd.read_csv(filename, header=None, delimiter="\t", error_bad_lines=False)
    events = df[[2]].squeeze()
    
    # -------------- Lights out --------------------
    LIGHTS_OUT = df[events.str.contains('|'.join(words1), case=False, na=False)]
    lights_out = str(LIGHTS_OUT[[1]].squeeze())
    lights_out = re.search(r'\b\d?\d:\d\d\b:\d\d\b', lights_out)
    
    if lights_out is None:
        df_out.loc[i]['LIGHTS_OFF'] = '.'
    else:
        df_out.loc[i]['LIGHTS_OFF'] = lights_out.group(0)
    
    # ---------------- Lights on -----------------------
    LIGHTS_ON = df[events.str.contains('|'.join(words2), case=False, na=False)]
    lights_on = str(LIGHTS_ON[[1]].squeeze())
    lights_on = re.search(r'\b\d?\d:\d\d\b:\d\d\b', lights_on)
    
    if lights_on is None:
        df_out.loc[i]['LIGHTS_ON'] = '.'
    else:
        df_out.loc[i]['LIGHTS_ON'] = lights_on.group(0)
    del df

df_out.to_csv('/Users/sq566/Desktop/wsc/logs/times.txt', sep ='\t', index=None)
