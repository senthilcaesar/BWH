function ChannelNumbers = getChannelNumbersFromHarmonizedLabels(HChannelSpreadsheet,Mrange,Labels,Transducers,Nmax)

%called from getEDFchannelnames()
%update to use settings.HChannelS struct
%no longer a need to import the whole spreadsheet

        ChannelNumbers = nan(length(Mrange),24);
        disp('Getting Harmonized Channel Label Info');
        HarmonizedChannelT1=readtable(HChannelSpreadsheet);
        HarmonizedChannelT = array2table(string(HarmonizedChannelT1{:,:}));
        HarmonizedChannelT.Properties.VariableNames = HarmonizedChannelT1.Properties.VariableNames;
        
        LabelsS=string(Labels);
        TransducersS=string(Transducers);
        LabelsTransducersS=strcat(LabelsS,'|',TransducersS);
        SigList = HarmonizedChannelT.Properties.VariableNames;
        for n=1:length(Mrange)
            if mod(n,20)==0
                disp(['n=' num2str(n)]);
            end
            temp = getChannelNumbersHarmonizedLabels([],HChannelSpreadsheet,LabelsS(n,:),TransducersS(n,:),LabelsTransducersS(n,:),HarmonizedChannelT{:,:},SigList);
            temp((Nmax+1):end)=[];
            ChannelNumbers(n,1:length(temp)) = temp;
        end

        ChannelNumbers = ChannelNumbers(:,1:Nmax);
        
        
        