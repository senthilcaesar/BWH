# Links
[Annotations table](#annotations)

[Signals table](#signals)

# Unmapped channel/annotation terms

_This page is for the internal use of the NSRR._

Ambiguous and currently unmapped terms encountered can be entered below. 

Upon review, the goal will be to choose one of three options for a given term:

 - ignore it

 - asign it to an existing mapping 
 
 - create a new mapping (i.e. a new primary annotation/channel concept that will be represented in the harmonized data)

If not ignoring it, there is a second choice to be made:

- assign it to the canonical/common files (`harm.annots` and `harm.canonical.sigs`)

- alternatively, assign it to a cohort specific mapping (`studies/{name}/harm.annots` or `studies/{study}/harm.sigs`)

If the term is idiosyncratic, conflicts with a canonical form, or otherwise confusing, then add to the study-specific file (or create that file).  When running NAP, these studys-specific mappings can over-ride the canonical mappings.

---


## Annotations

|   Term      | Study        | Example file    |    Notes                                            | Resolution     |
|   ----------| -------------| ----------------|    -------------------------------------------------| ---------------|
|   AE_Tech_Hy0% |  HCHS|   hchs-sol-15445599, hchs-sol-00302529|   Need to check if AE_Tech_Hy is the real hypopneas|    |
| AE_Tech_Hy1% |	HCHS|	hchs-sol-15445599, hchs-sol-00302530|	|    |
| AE_Tech_Hy3% |	HCHS|	hchs-sol-15445599, hchs-sol-00302531|	|    |
| AE_Tech_Hy4% |	HCHS|	hchs-sol-15445599, hchs-sol-00302532|	|    |
| Hypopnea-55-2-0 |	HCHS|	hchs-sol-15445599, hchs-sol-60285261|	|    |
| Hypopnea-55-2-1 |	HCHS|	hchs-sol-15445599, hchs-sol-60285261|	|    |
| Hypopnea-56-2-0 |	HCHS|	hchs-sol-15445599, hchs-sol-60285261|	|    |
| Hypopnea-3506-3-0 |	HCHS|	hchs-sol-38937437|	|    |
| Hypopnea-3507-3-0 |	HCHS|	hchs-sol-38937437|	|    |
| Apnea-55-1-0 |	HCHS|	hchs-sol-15445599,hchs-sol-38937437 |   |    |
| Apnea-55-1-1 |	HCHS|	hchs-sol-15445599, hchs-sol-38937437|	|    |
| Apnea-56-1-0 |	HCHS|	hchs-sol-15445599, hchs-sol-38937437|	|    |


## Signals

|   Term      | Study        | Example file    |    Notes                                            | Resolution     |
|   ----------| -------------| ----------------|    -------------------------------------------------| ---------------|
| A-Flow |  HCHS|hchs-sol-99999935  |   Nasal cannula signal in HCHS is marked as A-Flow. Need study specific case here.  |    |
| A-Effort/FVP |	HCHS|	|Ares Unicorder Effort-central venous pressure change	|    |
| Airflow | MrOS (visit 1, 2), SHHS(1,2), SOF, ABC (baseline,month09, month18),CCSHS,CFS, Chat(baseline, follow-up, non-randomized), HAASA, Homepap(Lab-full, lab-split) 	| 	|(1) Airflow is thermistor in MrOS, SHHS, SOF, ABC ,CCHS (thermistor/thermocouple), CFS, HAASA  (2) Chat airflow is Dymedix airflow with RERA-quality cable (3) Homepap airflow uspecified 	|     |
| Flow | MESA,Bestair (baseline,followup,non-randomized),Chat (baseline,followup,non-randomized), Heartbeat(baseline,followup),Homepap(home,lab-full,lab-split),NuMoM2b (visit 1,visit 3)| 	|flow signal may be the cannula pressure for all studies in MESA, unspecified flow in rest of datasets |     |
| Pres | MESA	| 	|Pres signal also present in MESA in addition to Flow, looks similar to flow	|     |
| Flow_EG|	Heartbeat(followup)|	|	|    |
| FLOW2|	Homepap (lab-full,lab-split)|	|	|    |
| FLOW5|	Homepap (lab-full,lab-split)|	|	|    |
| Flow_DC7|	Chat (non-randomized)|	|	|    |
| Flow_Patient|	Bestair (baseline)|	|	|    |
| Flow_Patient_1  | Bestair (baseline)	|	|	|    |
| Flow_Patient_2 | Bestair (baseline)|  |	|    |
| Flow_Patient_3 | Bestair (baseline)|	|	|    |
| XFlow |Bestair (baseline,followup,non-randomized),Chat (baseline,followup,non-randomized),Heartbeat(baseline,followup),Homepap (lab-full,lab-split)|  | Looks like some derived flow (Derived from Embla XactTrace belt in Bestair, Heartbeat, Homepap)	|    |
| XFlow_EG|	Heartbeat(followup)|	|	|    |
| XFlow_PDS|	Homepap(home)|	|	|    |
| NPV_flow |Homepap (lab-full)|	    |	|    |
| NPV |Homepap (lab-full)|  |	|    |
| New Air,new_air,New_Air,NEWAIR,NEW_AIR,New_A/F |SHHS 1|  |	looks like ProTech thermistor|    |
| PNEUMFLOW |	Chat (baseline)|    |	|    |
| Pressure |	Homepap (lab-split)|	 | May be CPAP pressure since its a split night study   |    |
| PressureE |	Homepap (lab-split)|	 |    |    |
| PressureI |	Homepap (lab-split)|	 |    |    |
| CFlow |	Bestair (baseline)|	bestair-baseline-400335| CFlow may represent either Cannula flow or CPAP flow depending on the specific study. CPAP signal seems on for this particular example study. Need to check study specific files to come up with mapping.	|    |
| Flo,Flow_A10,Flow_CU |	|	|Was present in the previous version of harm.annots. Couldn't trace to a specific dataset.	|    |
| Nasal/Oral |	|	|Was present in the previous version of harm.annots. Couldn't trace to a specific dataset.	|    |
| Resp_Airflow,Resp_Flow,Resp_FLOW-Ref |	|	|Was present in the previous version of harm.annots. Couldn't trace to a specific dataset.	|    |

