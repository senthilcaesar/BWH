# Running the Respiratory Endotype pipeline in NAP

## Input Data Requirements/Format

The estimation of respiratory endotypes is automated and executed
using custom software (Phenotyping Using Polysomnography (PUPbeta);
MATLAB, Mathworks, Natick MA) within the NAP pipeline. To run the
respiratory endotype program successfully, the user must be aware of
the following requirements for the input data:

##  Matlab Software & Toolboxes

Running PUPbeta requires Matlab 2019b or higher version with the following toolboxes installed:

 - Signal Processing
 - Curve Fitting
 - Statistics and Machine Learning
 - Optimization Toolbox
 - Wavelet
 - System Identification
 - Control Systems
 - Parallel Computing


##  EDF and XML/annotation File Information

The PUPbeta is run using the following command:

```
matlab -nodisplay \
       -r "StartHere ${resp_output} ${id}/ ${edfname} ${xmlname}" \
       -sd ${NAP_DIR}"/respiratory_endotypes" \
       -logfile ${resp_output}/${id}/outputconvert.log 2>> resp.log
```

Note that the program requires subject id (`id`), full path for EDF (`$edfname`), 
full path for XML files (`$xmlname`), full path for output folder (`$resp_output`) 
and the NAP code directory (`$NAP_DIR`).  For example:

```
id="100206525-201909-psg-DA7020MM"

edfname="/data/nsrr/working/puptest/may4/100206525-201909-psg-DA7020MM-rs1-harm.edf"

xmlname="/data/nsrr/working/puptest/may4/harm.annot"

resp_output="/data/nsrr/working/puptest/may4/nap/"

NAP_DIR='/PHShome/rma56/nsrr/nap'
```

## Signals and Annotations

### Nasal pressure

A nasal pressure signal (or flow from nasal cannula) must be present
in the EDF, since nasal pressure is used for determining the
ventilation and ventilatory drive. Without nasal pressure signal
present, respiratory endotype analysis should not be run.

!!!hint "Ideally the signal should be named as `csCAN` or `nas_pres`.

Care must be taken to map/harmonize only cannula signal to `csCAN` or
`nas_pres`. Thermistor or unknown airflow signal (i.e., airflow
channel without known transducer type or airflow derived from
thermistor) should not be mapped to `csCAN` or `nas_pres`.

### EEG 

Study should contain at least one EEG signal and cortical arousals
scored properly (arousals must be marked for the entire duration of
cortical arousals).

If the study does not contain at least one EEG signal and, or if
cortical arousals are not scored accurately (e.g., using a hot key to
denote arousal starting time with no duration information or using a
3s duration for all arousals), arousal related endotypic traits
(ArThresh, VRA) will not be estimated.  This could also lead to
overestimation of loop gain and other ventilatory parameters.

### SpO2

Another signal that must be present is the SpO2. This is used to
calculate respiratory event related desaturation, hypoxic burden
etc. within the software.
 
### Event annotations

Sleep and respiratory events must be provided in a single annotation
file in either XML, Luna's `.annot` or `.eannot` formats before
running PUPbeta.
 
NAP by default uses the `.annot` text format (described here).  The
tab-delimited columns should be

```
class     instance    channel     start stop  meta
```

Time stamps for `start` and `stop` shoudl be in either `HH:MM:SS.FFF`
or `HH:MM:SS` format.
 
#### Staging

For the sleep stages, coding required is either `N1`, `N2`, `N3`, `R`
and `W`, or `?` for unknown.  Alternate labels also accepted are
`wake`, `NREM1`, `NREM2`, `NMREM3`, `NREM4` and `REM`.

#### Respiratory events

For respiratory events the coding required: `arousal`, `apneacentral`,
`apneamixed`, `apneaobstructive`, and `hypopnea`.
 
In case of clinical studies where arousal/respiratory event scoring
are unreliable, PUPbeta can be used to auto score arousals and
respiratory events and use them for calculation of endotypes. Use the
below settings in `StartHere.m`.

```
settings.useWSanalysisToReplaceAr = 1;      % 0 = use original scoring,  
                                            % 1 = auto arousal scoring with best EEG

settings.AutoScoreRespEvents = 1;           % autoscoring of apneas and hypopneas

settings.UseAutoScoredRespEventsForLG = 1;  % 1= use autoscored respiratory events in endotyping
```
 
For running PUPbeta without autoscoring, one needs to turn the above
settings off (use 0) in `StartHere.m`.

### Body position

When called from within the NSRR NAP pipeline, we are not currently
considering the body position while calculating the respiratory
endotypes. The software assumes supine position for the entire study.

 