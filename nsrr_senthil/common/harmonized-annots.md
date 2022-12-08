# Harmonized annotations

These tables reflect the definitions in
```
nsrr/common/resources/harm.annots
```

which is the underlying file used by NAP to create the _harmonized_ annotations
(and should be kept in sync with the former).

Note the use of `class/instance` annotation naming (i.e. columns 1 and 2 of the .annot file). i.e.
this means that Luna masks will by default match all arousal types if just `arousal` is specified.

Cohort-specific definitions can over-ride these definitions by adding
additional definitions in a file (same format as `harm.annots`) and
passing the location to NAP via the `NAP_HARM_USER_ANNOTS`
configuration flag.

## Arousals

| Annotation                | Definition                              |
| ------------------------- | --------------------------------------- |
| arousal                   | EEG Arousal (generic; no specific subtype) |
| arousal/spontaneous       | Spontaneous EEG Arousal                     |                        
| arousal/respiratory       | Respiratory related EEG Arousal             |
| arousal/lm                | Limb movement related EEG arousal           |
| arousal/plm               | Periodic limb movement related EEG arousal  |
| arousal/external          | External EEG arousal |
| arousal/cheshire          | Arousal resulting from Chin EMG |


## Respiratory Events

Primary apnea and hypopnea

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| apnea                      | Apnea (unspecified type)                |
| apnea/obstructive          | Obstructive Apnea                       |
| apnea/central              | Central Apnea                           |
| apnea/mixed                | Mixed Apnea                             |
| hypopnea                   | Hypopnea (unspecified obstructive or central)|
| hypopnea/obstructive       | Obstructive Hypopnea                    |
| hypopnea/central           | Central Hypopnea                        |

Secondary:

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| periodic\_breathing        | Periodic breathing |
| cheynestokes\_breathing    | Cheyne Stokes breathing                 |
| respiratory\_paradox       | Respiratory (abdominal-thorax) paradox  |
| snoring                    | Snoring                                 |


## Oxygen Saturation Events 

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| desat                      | SpO2 desaturation                       |



## Staging (epoch-level annotations)

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| N1                      | Stage 1 sleep                           |
| N2                      | Stage 2 sleep                           |
| N3                      | Stage 3/4 sleep                         |
| R                       | REM sleep                               |
| W                       | Wake                                    |
| ?                       | Unknown stage                           |
| U                       | Unscored (effectively treated same as ?) |
| M                       | Body Movement                           |
| L                       | Lights On Epoch                         |
| lights\_on              | Lights On                               |
| lights\_off             | Lights Off                              |


## Limb Movements

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| movement                   | Generic movement (nb. epoch-level staging term is `M`) |
| | |
| PLM                        | Periodic Limb Movement |
| PLM/left                   | Periodic Limb Movement (Left)           |
| PLM/right                  | Periodic Limb Movement (Right)          |
| PLM/both                   | Periodic Limb Movement (Both limbs)     |
| | | 
| LM                         | Limb Movement (unspecified limb)        |
| LM/left                    | Limb Movement (left)                    |
| LM/right                   | Limb Movement (right)                   |
| LM/both                    | Limb Movement (both legs)               |

Note: previously `lml_arousal`, `lmr_arousal` and `lmb_arousal` were listed,
although these should be given by the `arousal/lm` above.


## Artifact

It is expected that artifact annotations will typically have an associated channel specifier.

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| artifact                   | Unspecified signal artifact             |
| artifact/respiratory       | Respiratory artifact                    |
| artifact/TcCO2             | Transcutaneous CO2 artifact             |
| artifact/SpO2              | SpO2 artifact                           |
| artifact/EtCO2             | End tidal CO2 artifact                  |
| artifact/blood\_pressure   | Blood pressure artifact  |
| artifact/proximal\_pH      | Proximal pH artifact |
| artifact/distal\_pH        | Distal pH artifact |
| artifact/pH	             | Unspecified pH artifact |
| artifact/body\_temperature | Body temperature artifact |
			


## Arrhythmias

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| arrhythmia/bradycardia                | Bradycardia                             |
| arrhythmia/tachycardia                | Tachycardia (unspecified)               |
| arrhythmia/narrow\_complex\_tachycardia | Narrow Complex Tachycardia              |


## Body position

Primary:

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| position/supine            | Supine position |
| position/nonsupine         | Non-supine position |

Secondary terms:

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| position/left              | Left side position | 
| position/right             | Right side position |
| position/upright           | Upright position |
| position/prone             | Prone position |
| position/unknown           | Explicitly indicating body position is unknown |


## Misc

| Annotation                 | Definition                              |
| -------------------------- | --------------------------------------- |
| notes | Technician notes |



