import pandas as pd
import re
import time
'''
This code read the annotation csv files and searches for
LIGHTS OFF and LIGHTS ON time parameter and save it was tab delimited text file.
Be aware of few discrepancies listed below.

Annotations discrepancy ( Apples )
-----------------------------------	
Lights OFF terms  | 	 Lights ON terms
-----------------------------------
ights out	      |    ?ights on
lgiths out	      |    liights on
loghts out	      |    ligts on
lighst out	      |    L-on
L-out	          |    lon
L-off	          |    l on
loff              |
l off             |	

'''
subjects = '/Users/sq566/Downloads/nsrr_studies/apples/data/subjects.txt'
with open(subjects) as f:
    case_list = f.read().splitlines()     
df_out = pd.DataFrame(index=range(0,len(case_list)), columns=['ID','LIGHTS_OFF'
                                                              ,'LIGHTS_ON'])

for i, ID in enumerate(case_list):

    df_out.loc[i]['ID'] = ID
    annot_file = f'/Users/sq566/Downloads/nsrr_studies/apples/data/csv/{ID}.csv'
    df = pd.read_csv(annot_file, encoding = "ISO-8859-1")
    df['Event type'] = df['Event type'].str.lower()
    df['Event type'] = df['Event type'].str.replace('supine','',flags=re.I)
    df['Time'] = pd.to_datetime(df['Time'], errors='coerce')
    
    # 1st filter
    searchfor1 = ['light', 'lights', 'loff', 'lon', 'lout', 'l off', 'l on', 'l out']
    words = [rf'\b{string}\b' for string in searchfor1]
    lights = df[df['Event type'].str.contains('|'.join(words), case=False, na=False)]
    print(ID)
    if not lights.empty:
        print(lights)
    time.sleep(1)
    
    # Lights off filter
    searchfor2 = ['out', 'off', 'ou', 'put']
    LIGHTS_OFF = lights[lights['Event type'].str.contains('|'.join(searchfor2))]
    if not LIGHTS_OFF.empty:
        lights_off = str(LIGHTS_OFF.Time.iloc[0])
        LOFF = re.search(r'\b\d?\d:\d\d\b:\d\d\b', lights_off)
        #print('found Lights off', LOFF.group(0))
        df_out.loc[i]['LIGHTS_OFF'] = LOFF.group(0)
    else:
      print('did not find lights off ', ID)
      df_out.loc[i]['LIGHTS_OFF'] = '.'
    
    # Lights on filter
    searchfor3 = ['on', 'in', 'oh']
    LIGHTS_ON = lights[lights['Event type'].str.contains('|'.join(searchfor3))]
    
    if not LIGHTS_ON.empty:
        lights_on = str(LIGHTS_ON.Time.iloc[-1])
        LON = re.search(r'\b\d?\d:\d\d\b:\d\d\b', lights_on)
        #print('found Lights on', LON.group(0))
        df_out.loc[i]['LIGHTS_ON'] = LON.group(0)
    else:
      print('did not find lights on ', ID)
      df_out.loc[i]['LIGHTS_ON'] = '.'
    
    print("\n")
    print("--------------------------------------------------------------------------")
    del df
      
df_out.to_csv('/Users/sq566/Downloads/nsrr_studies/apples/data/times.txt', sep ='\t', index=None) 
