% make tables for regression
%
% all the same number of rows
% matched on breath ID, which is subj and breath time
%
% PtData PtDataPnasal

%Breath ID
% ((Subj * 1E7) + BB_time) * 1000

BreathID = (PtData_flow.PT * 1E8) + (PtData_flow.BB_time*1E3);
length(unique(BreathID))
Atableflow = array2table([BreathID,Amatrix2_flow]);

BreathID_pnasal = (PtData_matched.PT * 1E8) + (PtData_matched.BB_time*1000);
length(unique(BreathID_pnasal))
pnasaltable = array2table([BreathID_pnasal,Amatrix2_pnasal_matched]);

nnz(ismember(BreathID, BreathID_pnasal))

combined_table = outerjoin(Atableflow, pnasaltable, 'key', 'Var1');


bbtimes_flow = PtData(PtData.PT==3,16);

bbtimes_flow = PtData_matched(PtData_matched.PT==3,16);
bbtimes_pnasal = PtData_pnasal_matched(PtData_pnasal_matched.PT==3,16);

nnz(bbtimes_flow.BB_time - bbtimes_pnasal.BB_time)

PtData(PtData.PT==3,:)


round(PtData_flow_pt.BB_time*1e3) - round(PtData_pnasal_pt.BB_time*1e3)

round(PtData_flow_pt.BB_time*1e3) - round(PtData{PtData.PT==3,16}*1e3)

length(round(PtData{PtData.PT==3,16}*1e3))


PT3 = PtData(PtData.PT==3,:);
PT3_bbindex = (PT3.PT * 1E8) + (PT3.BB_time*1000);
PT3 = [array2table(PT3_bbindex),PT3];

PT3_f = PtData_flow_pt;
PT3_bbindex = (PT3_f.PT * 1E8) + (PT3_f.BB_time*1000);
PT3_f = [array2table(PT3_bbindex),PT3_f];

PT3_p = PtData_pnasal_pt;
PT3_bbindex = (PT3_p.PT * 1E8) + (PT3_p.BB_time*1000);
PT3_p = [array2table(PT3_bbindex),PT3_p];

PT3_combined=outerjoin(PT3, PT3_f, 'key', 'PT3_bbindex');


Amatrix2_flow_pt
Amatrix2_pnasal_pt





%