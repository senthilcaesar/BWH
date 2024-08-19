# NSRR harmonized EDF signal

| Channel       | Reference            | Description                                     |
|---------------|----------------------|-------------------------------------------------|
| F3_M2         | F3 - M2              | A2/M2 (right mastoid EEG)                       |
| F3_LM         | F3 - (M1+M2)/2       |                                                 |
| F4_M1         | F4 - M1              | A1/M1 (left mastoid EEG)                        |
| F4_LM         | F4 - (M1+M2)/2       |                                                 |
| Fz_LM         | Fz - (M1+M2)/2       |                                                 |
| C3_M2         | C3 - M2              | C3 (left central EEG)                           |
| C3_LM         | C3 - (M1+M2)/2       |                                                 |
| C4_M1         | C4 - M1              | C4 (right central EEG)                          |
| C4_LM         | C4 - (M1+M2)/2       |                                                 |
| Cz_LM         | Cz - (M1+M2)/2       |                                                 |
| O1_M2         | O1 - M2              |                                                 |
| O1_LM         | O1 - (M1+M2)/2       |                                                 |
| O2_M1         | O2 - M1              |                                                 |
| O2_LM         | O2 - (M1+M2)/2       |                                                 |
| Oz_LM         | Oz - (M1+M2)/2       |                                                 |
| Pz_LM         | Pz - (M1+M2)/2       |                                                 |
| ECG           | ECG1 - ECG2          |                                                 |
| HR            |                      |                                                 |
| EMG           | EMG1 - EMG2          | EMG1 (Center chin emg), EMG2 (left submentalis emg) |
| LOC           | E1 - M2              | left eog                                        |
| ROC           | E2 - M2              | right eog                                       |
| thorax        |                      |                                                 |
| abdomen       |                      |                                                 |
| airflow       |                      | thermistor or a nasal pressure transducer       |
| thermistor    |                      |                                                 |
| nasal_pres    |                      |                                                 |
| SpO2          |                      | oximetry                                        |
| pleth         |                      | oximetry                                        |
| pulse         |                      |                                                 |
| LAT           | LLEG1 - LLEG2        |                                                 |
| RAT           | RLEG1 - RLEG2        |                                                 |
| cap           |                      | capnogram                                       |
| EtCO2         |                      | End-tidal CO2 measure derived from capnography  |
| TcCO2         |                      | Transcutaneous CO2                              |
| snore         |                      |                                                 |
| pos           |                      | body position                                   |
| cpap_flow     |                      | CPAP                                            |
| cpap_pressure |                      | CPAP                                            |


# NSRR harmonized annotation

| Harmonized class label | Description                                      |
|------------------------|--------------------------------------------------|
| N1                     |                                                  |
| N2                     |                                                  |
| N3                     |                                                  |
| W                      | Wake                                             |
| R                      | REM                                              |
| artifact               | SpO2 artifact, EtCO2 artifact, TcCO2 artifact, proximal pH artifact, Blood pressure artifact, Body_temperature_artifact |
| desat                  | SpO2 desaturation, Desaturation                  |
| arousal                | spontaneous arousal, respiratory arousal, arousal_lm, arousal_plm, Arousal_ASDA                |
| arrhythmia             | bradycardia, tachycardia                         |
| LM                     | Leg Movement, Limb Movement                      |
| PLM                    | Periodic Leg Movements                           |
| hypopnea               | obstructive, central                             |
| apnea                  | obstructive, central, mixed                      |
| biocal                 | eyes_blink, eyes_closed, eyes_open, inhale, exhale, flow_deep, flow_hold, flow_nasal, flow_oral, look_down, look_left, look_right, look_up, teeth_grind, cough, grit, hold_breath, mueller_maneuver, valsalva_maneuver, left_flex, right_flex |
| snoring                |                                                  |
| lights_off             |                                                  |
| lights_on              |                                                  |
| movement               |                                                  |
| pos2_supine            | supine                                           |
| pos2_nonsupine         | prone, left, right, upright                      |




