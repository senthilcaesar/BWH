function Pes_clean = PesRemoveEKG(Pes,EKG,dt,ECG_peak_i,templateLR,breakaftertemplate,Nbeatspersegment,contaminationmagnitudethreshold)
%load('J:\PEOPLE\FACULTY\SANDS\Phenotyping from Pes and Flow\SpikeFiles\1341.mat', 'Pes', 'EKG');
%Pes.values(end:end+1)=Pes.values(end);
%Pesclean = PesRemoveEKG(Pes.values,EKG.values,Pes.interval)
%addpath('E:\Work\MatlabFunctionDatabase')
%templateLR=[0.33, 1.67];
polyorder=3;
Time = (0:dt:(dt*(length(Pes)-1)))';
Niterations=3; %3
%Nbeatspersegment = 15;

if ~exist('ECG_peak_i')
    ECG_peak_i = EKGpeakdetection(EKG,Time,dt,0,0);
end
        medianRR = median(diff(ECG_peak_i))*dt;
        
        medianRRi = round(medianRR/dt);
    lefttemplatetime = templateLR(1)*medianRR;
    righttemplatetime = templateLR(2)*medianRR;
    leftdi1=round(1/dt*lefttemplatetime); rightdi1=round(1/dt*righttemplatetime);
    F1 = (leftdi1+rightdi1)/medianRRi;
    

   
[~,template,templateIQR] = crosscontaminatedPes(Pes,ECG_peak_i,leftdi1,rightdi1,1,polyorder,1); 
figure(89)
plot(template./nanstd(Pes)) 
hold on 
plot(templateIQR'./nanstd(Pes),'k:') 
w = hann(length(template));
plot(template(:).*w(:)./nanstd(Pes))
overallaveragecontamination = max(abs(template))/nanstd(Pes)
hold off
if breakaftertemplate
    return
end

% (contaminated,index,leftdi,rightdi,templateonly,polyorder,disppercentcomplete)


% clear temp contaminationmagnitude contaminationmagnitudepost
% processPes=1;
% loadprocessedPesifexists = 0; saveprocessedPes = 0; foundamatch=0;
% if loadprocessedPesifexists
%         foundamatch=0;
%         temp='Pes_clean';
%         for j=1:length(w)
%             foundamatch=strcmp(w(j).name,temp);
%             if foundamatch
%                 eval([temp '=filehandle.' temp ';']);
%                 break
%             end
%         end
%     if foundamatch
%         [~,template] = crosscontaminatedPes(Pes,ECG_peak_i,leftdi,rightdi,1);
%         contaminationmagnitude=max(abs(template));
%         [~,template] = crosscontaminatedPes(Pes_clean,ECG_peak_i,leftdi,rightdi,1);
%         contaminationmagnitudepost=max(abs(template));
%     end
% end
%%
% if (~foundamatch||loadprocessedPesifexists==0)&&processPes
try clf(111), catch me, end

Pes_clean=Pes;
plotalready=0;
    dividenightinto = floor(length(ECG_peak_i)/Nbeatspersegment);
    
for n=1:dividenightinto
    %%
    beatleft = Nbeatspersegment*(n-1)+1; beatright = beatleft + Nbeatspersegment + 1; %last 2 beats overlap
    
    if beatright>length(ECG_peak_i), beatright=length(ECG_peak_i); end
        ileft = ECG_peak_i(beatleft);
        iright = ECG_peak_i(beatright)+1;
    
    if n==1, ileft=1; end
    if n==dividenightinto, iright=length(EKG); end
    
    ECG_peak_i_segment = ECG_peak_i(beatleft:beatright) - ileft + 1;
    
    Pes_clean_segment=Pes_clean(ileft:iright);
    try
    leftdi=leftdi1; rightdi=rightdi1;
    
    for repeat=1:Niterations %repeats template subtraction with new template
        % EKG artifact removal
        [~,template] = crosscontaminatedPes(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,1,polyorder,0);
        figure(9)
        plot(template);
        if repeat==1
        contaminationmagnitude=max(abs(template))/nanstd(Pes_clean_segment)
        end
        
        %disp(['Decontaminating ' num2str(sum(contaminationmagnitude>contaminationmagnitudethreshold)) ' Pes signal']);
        if contaminationmagnitude>contaminationmagnitudethreshold
            Pes_clean_segment=crosscontaminatedPes(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,0,polyorder,0);
        end
        
        % after
        %[~,template] = crosscontaminatedPes(Pes_clean_segment,ECG_peak_i_segment,leftdi,rightdi,1);
        %contaminationmagnitudepost(i)=max(abs(template));
        if repeat<Niterations
            leftdi=round(leftdi*1.1); rightdi=round(rightdi*1.1);
            leftdi=leftdi+leftdi; rightdi=rightdi-leftdi;
        end
    end
    catch me
        disp(me.message)
    end
    Pes_clean(ileft:iright)=Pes_clean_segment;
    if 1
        %if ~plotalready
        figure(111)
        if ~plotalready
            axx(1)=subplot(2,1,1); plot(Time,EKG);
                axx(1).FontSmoothing = 'off';
            axx(2)=subplot(2,1,2); plot(Time,Pes,'k'); hold('on');
                axx(2).FontSmoothing = 'off';
            plotalready=1;
            linkaxes(axx,'x');
            xlim([Time(ileft)-5 Time(iright)+5]);
        end
        axx(2)=subplot(2,1,2); plot(Time(ileft:iright),Pes_clean(ileft:iright),'r');
        xlim([Time(ileft)-5 Time(iright)+5]);
    n
    end
end
%     if saveprocessedPes
%         currentdir = cd;
%         savestring = ['Pes_clean'];
%         temp = '-append';
%         cd(directory);
%         eval(['save ' filename ' ' savestring ' ' temp]);
%         cd(currentdir);
%     end
% else %no decontamination needed, just copy Pes into Pes_clean
%     Pes_clean=Pes;
% end

%%
if 1
    try clf(111), catch me, end
    figure(111)
    axx(1)=subplot(2,1,1); plot(Time,EKG);
    axx(2)=subplot(2,1,2); plot(Time,Pes,'r',Time,Pes_clean,'b');
    linkaxes(axx,'x');
end

