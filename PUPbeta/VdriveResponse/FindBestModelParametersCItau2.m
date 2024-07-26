function [Error,Vchem,Varousal,LoopGainParameters,BestParameters,BestSSres,FitQuality,i_1,i_end,lowerSEM,upperSEM,CI_parameters] = ...
    FindBestModelParametersCItau2(Data,VraOn,Veupnea,polyfitorder)
% Written by Philip I Terrill and Scott A Sands
% This MATLAB function accompanies the manuscript entitled:
% "Quantifying the Ventilatory Control Contribution to Obstructive Sleep Apnea"
% **************************************************************************
% **************************************************************************
%
% This function is used to identify the set of parameters (LG_0, Tau,
% delay, Vdr[1]) that results in the best fit between the first order
% ventilation model and the observed ventilation. The MATLAB function
% fmincon.m is used to search for the optimal parameter set.
% >>>>>>>>>>>>>>>>>>>>>>>>Input parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% ***Data: This matrix contains information about the patients observed
% ventilation
%             - Column 1: Breath-by-Breath Minute Ventilation -- Must already be "mean subtracted" i.e. eupnea=0.
%             - Column 2: Obstructive respiratory event series. Elements=1 for an un-obstructed breath; Elements=0 for an obstructed breath.
%             - Column 3: Arousal Scoring series: defines whether there is an arousal present in each breath (arousal present, AR[n]=1; No arousal present, AR[n]=0);
%             - Column 4: The time that each breath in the series occurs.
%             - Column 5: Scored Central Apnea array series: defines whether each breath is scored as a central apnea. elements=1 for breaths in a scored central apnea; Elements=0 for breaths not classified as a central apnea.
%             - Column 6: T_tot series (currently unused)
%             - Column 7: T_tot series for the previous breath.
% ***VraOn: This binary variable allows the inclusion or exclusion of the
% ventilatory response to arousal parameter.
% ***Veupnea: This is the eupneic ventilation level. This is typically normalised to 1.
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
% >>>>>>>>>>>>>>>>>>>>>>>>Output parameters<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% - Error: 1*N Array. This is the breath-by-breath different between the
% model predicted ventilation and the observed ventilation
% - Vchem: 1*N Array. The ventilatory drive due to chemical drive
% - Varousal: 1*N Array. The ventilatory drive due to arousal
% - LoopGainParameters: 1*8 Array of the derived loop gain parameters:
%           - Element 1: Steady State Loop Gain
%           - Element 2: Time Constant
%           - Element 3: Delay
%           - Element 4: Loop gain at the natural frequency
%           - Element 5: The natural cycling period
%           - Element 6: Loop gain at 1 cycle per minute
%           - Element 7: Loop gain at 2 cycles per minute
%           - Element 8: Additional ventilatory drive during breaths
%           occuring during an arrousal.
% - BestParametersForEachDelay: 1*5 Array of the fitted ventilatory
% control model parameters:
%           - Element 1: Steady State Loop Gain
%           - Element 2: Time Constant
%           - Element 3: Delay
%           - Element 4: Error in observed ventilattion to calculate Vchem_1
%           - Element 5: Additional ventilatory drive during breaths
%           occuring during an arrousal.
% - BestSSres: 1*1 Array. Mean Square Error (SSres) of the model fit to observed
% ventilation data.

% - FitQuality: 1*4 Array. Summary Measures of the quality of the fit between
% the model and observed ventilation
%           - Element 1: Mean Square Error (SSres) of the model fit to observed
%           ventilation data.
%           - Element 2: R_squared for the correlation between the observed
%           ventilation and model predicted ventilation.
% - i_1: The index of the first breath (n=1) in the ventilation series for which Vchem and Vdrive are estimated
% - i_end: The index of the last breath in the ventilation series for which
% the model was applied
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%
%**************************************************************************
%**************************************************************************
global settings
ploton = settings.plotfigureLGparameters;

if ~isfield(settings,'tau2ontau1')
    settings.tau2ontau1 = 0;
end

warning('off');
Minute_Ventilation=Data(:,1); % The Minute Ventilation array is a breath-by-breath minute ventilation series which is normalised
ObstructiveEvents=Data(:,2); % The ObstructiveEvents array defines whether each breath is obstructed. Penalize[n]=1 for an un-obstructed breath; Penalize[n]=0 for an obstructed breath.
TimeB=Data(:,4); % The time of each breath, from the first breath in the series.
Ttot_ser=Data(:,6); % T_tot Series
Ttot_median=median(Ttot_ser); % The median t-tot based on the T_tot series [diff(TimeB)].


%%
%**************************************************************************
% Set the constraints for the optimisation procedure:


%i_1=i_min; % The index of the first breath for which the model will be fit.
i_end=length(Minute_Ventilation);

%**************************************************************************
% Count the number of obstructive events within the analysis window. If
% there is not at least one obstructive event commencing within the window,
% pass all output parameters as NaN's
%**************************************************************************

% if  N_events<1
%     LoopGainParameters(1:8)=NaN;
%     FitQuality(1:3)=NaN;
%     Error=NaN;
%     Vchem=NaN;
%     Varousal=NaN;
%     BestSSres=NaN;
%     BestParameters=NaN;
%     i_1=NaN;
%     i_end=NaN;
%     lowerSEM=NaN;
%     upperSEM=NaN;
%     CI_parameters(1:5,1)=NaN;
%     return
% end
%**************************************************************************
%**************************************************************************
i_1=10; % The index of the first breath for which the model will be fit.

secresolution=1;
mindelaysec=3;
mindelayi=round(mindelaysec/secresolution);

delta_delay = ceil(1/(Ttot_median/secresolution));

maxdelaysec = 20;
maxdelayi = round(maxdelaysec/Ttot_median);

%if delta_delay>1, delta_delay=1; end
%delay_breaths=[delta_delay:delta_delay:7];
delay_breaths=(mindelayi*delta_delay):delta_delay:maxdelayi;

secresolutionactual=delta_delay*Ttot_median;

N_events=0;
for k=i_1:length(ObstructiveEvents)
    if ObstructiveEvents(k)==0&&ObstructiveEvents(k-1)==1
        N_events=N_events+1;
    end
end

%%

for dd=1:length(delay_breaths)
    %**********************************************************************
    % Generate the delayed VE:
    %**********************************************************************
    delay_secs=delay_breaths(dd)*Ttot_median;
    delayed_VE{dd}=zeros(length(Minute_Ventilation),1);
    for k=i_1:length(Minute_Ventilation)
        d_i=Ttot_ser(k-1);
        k_back=1;
        while delay_secs>d_i
            k_back=k_back+1;
            ser_range = k-k_back:k-1;
            if (any(ser_range==0)) % if we go too far
                d_i=delay_secs; % then set d_i
                k_back=k_back-1; % and k_back, to avaoid indexing problems
                break
            else
                d_i=sum(Ttot_ser(ser_range));
            end
        end
        delayed_VE{dd}(k)=Minute_Ventilation(k-k_back)+(1-(delay_secs-sum(Ttot_ser(k-k_back+1:k-1)))/Ttot_ser(k-k_back))*(Minute_Ventilation(k-k_back+1)-Minute_Ventilation(k-k_back));
        delayed_VE{dd}(1:i_1-1)=-1; %why this, not zero
        delayed_VE{dd}(1:i_1-1)=0; %why this, not zero
    end
end
%**********************************************************************
%**********************************************************************

%**************************************************************************
% In order to identify the optimal parameter set (of LG_0, tau, Vdr[1],
% VRA) the inbuilt MATLAB function fmincon.m is utilised with an interior
% point algorithm.
%
% To maximize the robustness of the fit (avoiding local minima), a set of
% five starting conditions (initial parameter values) were applied. The
% optimisation algorithm is applied for up to 400 iterations.
%
% The starting set that yielded the the lowest final SSres is identified. The
% final end parameter set for this optimal starting point is now used for
% up to 1000 further iterations to identify the final optimal paramter set:
%**************************************************************************
%%

fLG=12;
ftau=36;
fVRA=20;
%
% fLG=2;
% ftau=2;
% fVRA=20;
%fdont=5;

sLG=0.125;
stau=5;
sVRA=0.05;
%sdont=1/3;

%initial limits
mLG=2;
mtau=200;
mVRA=1; %fraction of Veupnea
%mdont=3;

%hard limits
maxLG1=3;
maxtau=360;
maxVRA=1; %fraction of Veupnea
mintau=3;

Parameters=0;
LG1test=sLG*fLG.^[0:20]; LG1test(LG1test>mLG)=[];
%LG1test=sLG;
tau1test=stau*ftau.^[0:20]; tau1test(tau1test>mtau)=[];
VRAtest=Veupnea*[sVRA*fVRA.^[0:20]]; VRAtest(VRAtest>mVRA)=[];
%VRAtest=Veupnea*sVRA;
if VraOn==0
    VRAtest=VRAtest*0.00000000001;
end
%delayonsqrttau1test=sdont*fdont.^[0:20]; delayonsqrttau1test(delayonsqrttau1test>mdont)=[];

SSres_LG1 = NaN*ones(length(LG1test),length(tau1test),length(VRAtest));
value = NaN*zeros(length(delay_breaths),1);
SSres_LG1_dd = NaN*ones(length(LG1test),length(tau1test),length(VRAtest),length(delay_breaths));

%delayplussqrttau = delay_breaths*Ttot_median*3;

dd_range = 1:length(delay_breaths);
for m=1:length(delay_breaths)
    dd = dd_range(m);
    %OPTIONS = optimset('Display','off','TolX',1e-4,'TolFun',1e-4,'MaxIter',5,'MaxFunEvals',5,'Algorithm','interior-point');
    clear SSres_LG1
    for i=1:length(LG1test)
        for j=1:length(tau1test)
            for k=1:length(VRAtest)
                if 1
                    %delay_z=delayonsqrttau1test(j)*delayplussqrttau(dd)/(delayonsqrttau1test(j)+1);
                    %                     delay_z=delayonsqrttau1test(j)*delayplussqrttau(dd)/(delayonsqrttau1test(j)+1);
                    %                     delay_zi=round((delay_z/secresolutionactual));
                    %                     if delay_zi<1, delay_zi=1; end
                    %                     if delay_zi>length(delay_breaths), delay_zi=length(delay_breaths); end
                    %                     delay_z=delay_zi*secresolutionactual;
                    %                     tau_z=(delayplussqrttau(dd)-delay_z)^2;
                    %                     if tau_z<Ttot_median/4, tau_z=Ttot_median/4; end
                    %                     dd1 = delay_zi;
                    %[final_parameters_LG1(i,j,k),SSres_LG1(i,j,k),~,~] = fmincon(@(Parameters) TheModelX1_knownLG1tau(Parameters,Data,i_1,delayed_VE{dd},polyfitorder,delay_secs,LG1test(i),tau1test(j),VRAtest(k)),Parameters,[],[],[],[],lower(3),upper(3),@(Parameters) model_nonlinear_constraints(),OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
                    ParametersFixed = [LG1test(i),tau1test(j),NaN,VRAtest(k)];
                    [SSres_LG1(i,j,k),~,~] = TheModelX1tau2(Parameters,Data,i_1,delayed_VE{dd},polyfitorder,ParametersFixed);
                end
            end
        end
    end
    [value(dd),I] = min(SSres_LG1(:));
    [I1,I2,I3] = ind2sub(size(SSres_LG1),I);
    if 1
        LG1test_(dd)=LG1test(I1); %temp(I1,1,I2,I3);
        tau1test_(dd)=tau1test(I2); %temp(I1,2,I2,I3); ................ up to here..
        %delayonsqrttau1test_(dd)=delayonsqrttau1test(I2); %temp(I1,2,I2,I3);
        VRAtest_(dd)=VRAtest(I3); %temp(I1,4,I2,I3);
        if 0
            delay_z=delayonsqrttau1test_(dd)*delayplussqrttau(dd)/(delayonsqrttau1test_(dd)+1);
            delay_zi=round((delay_z/secresolutionactual));
            if delay_zi<1, delay_zi=1; end
            if delay_zi>length(delay_breaths), delay_zi=length(delay_breaths); end
            delay_z=delay_zi*secresolutionactual;
            delaytest_(dd)=delay_z/(Ttot_median);
            tau1test_(dd)=(delayplussqrttau(dd)-delay_z)^2;
            if tau1test_(dd)<Ttot_median/4, tau1test_(dd)=Ttot_median/4; end
        end
    end
    %if ~excludethis
    SSres_LG1_dd(:,:,:,dd) = SSres_LG1;
    %end
end
[~,bestdd1]=min(value); bestdd11=bestdd1;

[minvalueall,I] = min(SSres_LG1_dd(:));
[besti,bestj,bestk,bestdd] = ind2sub(size(SSres_LG1_dd),I);
if ploton
    figure(21); clf(figure(21)); set(gcf,'color',[1 1 1]);
    loglog(delay_breaths*Ttot_median,value,'k.-'); box('off'); hold('on');
    loglog(LG1test,SSres_LG1_dd(:,bestj,bestk,bestdd1),'b.-');
    loglog(tau1test/60,squeeze(SSres_LG1_dd(besti,:,bestk,bestdd1)),'r.-');
    loglog(VRAtest+0.001,squeeze(SSres_LG1_dd(besti,bestj,:,bestdd1)),'g.-');
    %hold('off');
end

%all breaths rerun
delay_breaths2 = delay_breaths;
ddr=length(delay_breaths);
ddl=1;
dd_range=ddl:ddr;

I_remove = 1+0*delay_breaths; I_remove(dd_range)=0;

LG1test_(I_remove==1) = NaN;
%delayonsqrttau1test_(I_remove==1) = NaN;
VRAtest_(I_remove==1) = NaN;
%Error_(I_remove==1) = NaN;
%delaytest_(I_remove==1) = NaN;
tau1test_(I_remove==1) = NaN;
value(I_remove==1) = NaN;
clear I1 I2 I3

%% Step 2
closeinby1=0.5;
%besterror=[];
closeinbyrate=0.5;

%%
%to do: move delay range around to meet best delay +/- 1.
closeinby=closeinby1; %current grid size
SSres_LG1_dd = NaN*ones(3,3,3,length(delay_breaths));
nn_thres=[99 99 5 2 1.5*ones(1,10)];
value_nn=[];
for nn=1:length(nn_thres) %progress loop, move grid or zoom
    
    %Should we zoom in (else: "walk" the grid)? (was the best result at the middle of last grid besti==2,bestj==2,bestj==3)
    %ok to zoom if term is at its own limit since we don't want to walk the grid further in that direction.
    if nn>1&&(besti==2||LG1test_last(bestdd11)<0.01||LG1test_last(bestdd11)>=maxLG1)&&(bestj==2||tau1test_last(bestdd11)<mintau||tau1test_last(bestdd11)>=maxtau)&&(bestk==2||VRAtest_last(bestdd11)<(0.01*Veupnea)||VRAtest_last(bestdd11)>=(maxVRA*Veupnea)) %%soft limits are effectively set here
        closeinby=closeinby*closeinbyrate;
    else
        %do nothing -- walk grid by current step size (closeinby stays the same)
    end
    
    LG1test_last=LG1test_; tau1test_last=tau1test_; VRAtest_last=VRAtest_; %for plotting
    
    for m=1:length(delay_breaths2) %set of possible delays
        dd = dd_range(m);
        
        % Should we skip this delay?
        % yes if results were very poor: more likely to skip in as loops go on [99 -> 1.5]
        % but never skip neighbors of best delay (dd>bestdd1+1||dd<bestdd1-1)
        if value(dd)>nn_thres(nn)*value(bestdd)&&(dd>bestdd1+1||dd<bestdd1-1)
            continue
        end
        
        %         if dd<bestdd1-nn_thres(nn)||dd>bestdd1+nn_thres(nn)
        %             continue
        %         end
        Parameters=0;
        clear SSres_LG1
        LG1test=LG1test_(dd)*[1/fLG 1 fLG].^closeinby; %[1/fLG 1 fLG].^closeinby is the latest grid for search
        %         if tau1test_(dd)>10*180
        %             tau1test_(dd)=10*180;
        %         end
        tau1test=tau1test_(dd)*[1/ftau 1 ftau].^closeinby; %start closing in later after finding better delays.
        tau1test(tau1test>maxtau)=maxtau;
        tau1test(tau1test<mintau)=mintau;
        %         if VRAtest_(dd)>Veupnea
        %             VRAtest_(dd)=Veupnea;
        %         end
        VRAtest=VRAtest_(dd)*[1/fVRA 1 fVRA].^closeinby;
        VRAtest(VRAtest>(maxVRA*Veupnea))=(maxVRA*Veupnea);
        
        %         if VraOn==0 %overwrite
        %             VRAtest=0;
        %         end
        for i=1:length(LG1test)
            for j=1:length(tau1test)
                for k=1:length(VRAtest)
                    if 1
                        %                         delay_z=delayonsqrttau1test(j)*delayplussqrttau(dd)/(delayonsqrttau1test(j)+1); %sec
                        %                         delay_zi=round((delay_z/secresolutionactual));
                        %                         if delay_zi<1, delay_zi=1; end
                        %                         if delay_zi>length(delay_breaths), delay_zi=length(delay_breaths); end
                        %                         delay_z=delay_zi*secresolutionactual;
                        %                         tau_z=(delayplussqrttau(dd)-delay_z)^2;
                        %                         if tau_z<Ttot_median/4, tau_z=Ttot_median/4; end
                        %                         dd1 = delay_zi;
                        %[final_parameters_LG1(i,j,k),SSres_LG1(i,j,k),~,~] = fmincon(@(Parameters) TheModelX1_knownLG1tau(Parameters,Data,i_1,delayed_VE{dd},polyfitorder,delay_secs,LG1test(i),tau1test(j),VRAtest(k)),Parameters,[],[],[],[],lower(3),upper(3),@(Parameters) model_nonlinear_constraints(),OPTIONS); %@(Parameters) model_nonlinear_constraints() is now [],[]
                        ParametersFixed = [LG1test(i),tau1test(j),NaN,VRAtest(k)];
                        [SSres_LG1(i,j,k),~,~] = TheModelX1tau2(Parameters,Data,i_1,delayed_VE{dd},polyfitorder,ParametersFixed);
                    end
                end
            end
        end
        [value(dd),I] = min(SSres_LG1(:));
        [I1,I2,I3] = ind2sub(size(SSres_LG1),I);
        if 1
            LG1test_(dd)=LG1test(I1); %temp(I1,1,I2,I3);
            tau1test_(dd)=tau1test(I2); %temp(I1,2,I2,I3);
            %delayonsqrttau1test_(dd)=delayonsqrttau1test(I2); %temp(I1,2,I2,I3);
            VRAtest_(dd)=VRAtest(I3); %temp(I1,4,I2,I3);
            %tau1test_(dd)=delayonsqrttau1test(I2); %temp(I1,2,I2,I3);
            if 0
                delay_z=delayonsqrttau1test_(dd)*delayplussqrttau(dd)/(delayonsqrttau1test_(dd)+1);
                delay_zi=round((delay_z/secresolutionactual));
                if delay_zi<1, delay_zi=1; end
                if delay_zi>length(delay_breaths), delay_zi=length(delay_breaths); end
                delay_z=delay_zi*secresolutionactual;
                delaytest_(dd)=delay_z/(Ttot_median);
                tau1test_(dd)=(delayplussqrttau(dd)-delay_z)^2;
                if tau1test_(dd)<Ttot_median/4, tau1test_(dd)=Ttot_median/4; end
            end
        end
        %if ~excludethis
        SSres_LG1_dd(:,:,:,dd) = SSres_LG1;
        %end
    end
    [temp,bestdd1]=min(value);
    %besterror=[besterror temp];
    [value_nn(nn),I] = min(SSres_LG1_dd(:)); %find minimum of 4D array by making it into 1D
    [besti,bestj,bestk,bestdd] = ind2sub(size(SSres_LG1_dd),I); %then lookup answer
    
    if ploton
        figure(21); set(gcf,'color',[1 1 1]);
        loglog(delay_breaths*Ttot_median,value,'k.-'); box('off'); hold('on');
        xvalues = LG1test_last(bestdd11)*[1/fLG 1 fLG].^closeinby;
        xvalues(xvalues>maxLG1)=maxLG1;
        loglog(xvalues,SSres_LG1_dd(:,bestj,bestk,bestdd1),'b.-');
        xvalues = tau1test_last(bestdd11)*[1/ftau 1 ftau].^closeinby/60;
        xvalues(xvalues>maxtau)=maxtau;
        xvalues(xvalues<mintau)=mintau;
        loglog(xvalues,squeeze(SSres_LG1_dd(besti,:,bestk,bestdd1)),'r.-');
        xvalues = VRAtest_last(bestdd11)*[1/fVRA 1 fVRA].^closeinby;
        xvalues(xvalues>(maxVRA*Veupnea))=(maxVRA*Veupnea);
        xvalues(xvalues==0)=0.001;
        loglog(xvalues,squeeze(SSres_LG1_dd(besti,bestj,:,bestdd1)),'g.-');
    end
    bestdd11=bestdd1;
    if nn>8
        if value_nn(nn)>value_nn(nn-4)/1.5 %
            break
        end
    end
end
%bestdd1 = delay_breaths2(bestdd1i);
%besterror

%%
if 0 %beware that you can't run a delay that hasn't been assessed above...
    delay_breaths3 = delay_breaths2;
elseif 0
    value2=value;
    value2(bestdd1)=NaN;
    [~,bestdd2]=min(value2);
    temp1 = sort([bestdd1 bestdd2]);
    ddr=temp1(2);
    ddl=temp1(1);
    delay_breaths3=delay_breaths(ddl:ddr);
elseif 0
    delay_breaths3 = delay_breaths(bestdd1);
    ddl=bestdd1;
    ddr=bestdd1;
else %keep best delay +/- 1 to go forward:
    ddr=bestdd1+1;
    ddl=bestdd1-1;
    if ddl<1;ddl=1; end
    if ddr>length(delay_breaths);ddr=length(delay_breaths); end
    delay_breaths3=delay_breaths(ddl:ddr);
end
dd_range=ddl:ddr;
I_remove = 1+0*delay_breaths; I_remove(dd_range)=0;
LG1test_(I_remove==1) = NaN;
tau1test_(I_remove==1) = NaN;
%delaytest_(I_remove==1) = NaN;
%delayonsqrttau1test_(I_remove==1) = NaN;
VRAtest_(I_remove==1) = NaN;
%value(I_remove==1) = NaN;
clear I1 I2 I3

% %% check/debug
% delay_z=(delayonsqrttau1test_.*delayplussqrttau./(delayonsqrttau1test_+1))/Ttot_median;
% delay_zB=(delayonsqrttau1test_.*delayplussqrttau./(delayonsqrttau1test_+1))/secresolutionactual

%%
maxiterN = 1000;
OPTIONS = optimset('Display','off','TolX',1e-8,'TolFun',1e-8,'MaxIter',maxiterN,'MaxFunEvals',maxiterN,'Algorithm','interior-point');
%%
Optimal_SSres = NaN*delay_breaths;
for m=1:length(delay_breaths3)
    dd = dd_range(m);
    start2 = [LG1test_(dd) tau1test_(dd) 0 VRAtest_(dd)];
    
    upper2 = start2*1.4;
    lower2 = start2/1.4;
    
    %simulation starting error range
    upper2(3)=start2(3)+Veupnea/2;
    lower2(3)=start2(3)-Veupnea/2;
   
    %handle boundaries:
        if upper2(4)>=(maxVRA*Veupnea)
            upper2(4)=(maxVRA*Veupnea);
            lower2(4)=(maxVRA*Veupnea)*0.9;
            start2(4)=(maxVRA*Veupnea);
        end
        
        if upper2(2)>=maxtau
            upper2(2)=maxtau;
            lower2(2)=maxtau*0.9;
            start2(2)=maxtau;
        end
        
        if upper2(1)>=maxLG1
            upper2(1)=maxLG1;
            lower2(1)=maxLG1*0.9;
            start2(1)=maxLG1;
        end
    
    ParametersFixed = NaN*start2;
    if VraOn==0
        start2(4)=[];
        lower2(4)=[];
        upper2(4)=[];
        ParametersFixed=[NaN NaN NaN 0];
    else
        ParametersFixed = NaN*start2;
    end
    
    Parameters = start2;
    % Use the end-point of the best starting condition as the start point for a further 1000 iterations:
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    if sum(isnan([upper2 lower2]))>0
        disp(num2str([upper2 lower2]));
    end
    
    %     dd1=round((delaytest_(dd)*Ttot_median)/secresolutionactual);
    %     if dd1<1, dd1=1; end
    %     if dd1>length(delayplussqrttau),dd1=length(delayplussqrttau); end
    %
    %     [Parameters dd1]
    [BestParametersForEachDelay(:,dd),Optimal_SSres(dd),~,~] = fmincon(@(Parameters) TheModelX1tau2(Parameters,Data,i_1,delayed_VE{dd},polyfitorder,ParametersFixed),Parameters,[],[],[],[],lower2,upper2,[],OPTIONS);
    
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    %**************************************************************************
    %**************************************************************************
end
[~, ind]=min(Optimal_SSres);
if ploton
    figure(21); set(gcf,'color',[1 1 1]);
    
    loglog(delay_breaths(ind)*Ttot_median,Optimal_SSres(ind),'k.','markersize',15); box('off'); hold('on');
    loglog(BestParametersForEachDelay(1,ind),Optimal_SSres(ind),'b.','markersize',15);
    loglog(BestParametersForEachDelay(2,ind)/60,Optimal_SSres(ind),'r.','markersize',15);
    if VraOn
        loglog(BestParametersForEachDelay(4,ind)+0.001,Optimal_SSres(ind),'g.','markersize',15);
    else
        loglog(0.001,Optimal_SSres(ind),'g.','markersize',15);
    end
    ylimtemp=get(gca,'ylim');
    ylim([Optimal_SSres(ind) min([Optimal_SSres(ind)*5,ylimtemp(2)])]);
    %loglog(delaytest_(ind)*Ttot_median,Optimal_SSres(ind),'ko','markersize',5); box('off'); hold('on');
    %loglog(BestParametersForEachDelay(2,ind)/60,Optimal_SSres(ind),'ro','markersize',5);
    hold('off');
end

%%

%ind-dd_range(2)
BestParameters=BestParametersForEachDelay(:,ind);
delay_secs=delay_breaths(ind)*Ttot_median; %%%%%%%%%%%%
if VraOn==0
    ParametersFixed=[NaN NaN NaN 0];
else
    ParametersFixed = NaN*BestParameters;
end
%dd1=round((delaytest_(ind)*Ttot_median)/secresolutionactual);
[BestSSres,Error,Rsq,Vchem,Varousal] = TheModelX1tau2(BestParameters,Data,i_1,delayed_VE{ind},polyfitorder,ParametersFixed);
FitQuality=[BestSSres Rsq delay_breaths(ind)]; % This array compiles the parameters describing the quality of the model fit to observed ventilation data
BestParameters(5)=delay_breaths(ind); % The optimal delay (breaths)

%Currently will perform Jacobian and 95%CI for best delay only:
noJacobian=0;
if ~noJacobian
    clear modelY_ modelY1 modelY2 f_ f1 f2
    %dX=[0.001 0.001 0.001 0.001 Ttot_median];
    dX=[0.1 1 0.1 0.1 Ttot_median];
    BestParametersDelay = BestParameters;
    BestParametersDelay(5) = delay_secs;
    for i=1:length(BestParametersDelay)
        parametersX1=BestParametersDelay;
        parametersX2=BestParametersDelay;
        parametersX1(i)=BestParametersDelay(i)+dX(i);
        parametersX2(i)=BestParametersDelay(i)-dX(i);
        [FF,Error_,~,VchemA,~,Penalized_Error] = TheModelX1tau2(BestParameters,Data,i_1,delayed_VE{ind},polyfitorder,NaN*BestParameters);
        modelY_=VchemA;
        [f1(i),~,~,VchemA,~] = TheModelX1delay(parametersX1,Data,i_1,polyfitorder);
        modelY1(:,i)=VchemA;
        [f2(i),~,~,VchemA,~] = TheModelX1delay(parametersX2,Data,i_1,polyfitorder);
        modelY2(:,i)=VchemA;
        
        if 0
            %Using minimum since large errors occur e.g. due to issues with 2delay<Tn<4delay bounds.
            JACOBIANA = 2*(modelY1(:,i)-modelY_)/dX(i);
            JACOBIANB = -2*(modelY2(:,i)-modelY_)/dX(i);
            temp = [(JACOBIANA) (JACOBIANB)];
            [~,tempI] = min([nanstd(JACOBIANA) nanstd(JACOBIANB)]);
            JACOBIAN(:,i) = temp(:,tempI);
        else
            %JACOBIAN(:,i) = -abs((modelY2(:,i)-modelY1(:,i))/dX(i));
            JACOBIAN(:,i) = ((modelY2(:,i)-modelY1(:,i))/dX(i));
        end
        RESIDUAL = Error_';
    end
    %     xdata=delayed_VE{ind}';
    %     [modelVdrive,delta] = nlpredci(@(BestParameters,xdata) TheModelX2(BestParameters,Data,i_1,xdata,polyfitorder,delay_secs),xdata,BestParameters,RESIDUAL,'Jacobian',JACOBIAN);
    if 0
        figure();
        plot(Data(:,4),Vchem); hold('on');
        plot(Data(:,4),modelY1(:,i),'color',[0.6 0.6 0.6]);
        plot(Data(:,4),modelY2(:,i),'color',[0.6 0.6 0.6]);
    end
    xdata=Data(:,1)';
    [modelVdrive,delta] = nlpredci(@(BestParametersDelay,xdata) TheModelX2delay(BestParametersDelay,Data,i_1,xdata,polyfitorder),xdata,BestParametersDelay,RESIDUAL,'Jacobian',JACOBIAN);
    
    if 0
        figure();
        plot(Data(:,4),Vchem); hold('on');
        plot(Data(:,4),modelVdrive+delta,'color',[0.6 0.6 0.6]);
        plot(Data(:,4),modelVdrive-delta,'color',[0.6 0.6 0.6]);
    end
    
    if size(modelVdrive,1)~=1 %to do -- work out why modelVdotMax is sometimes sideways...
        modelVdrive=modelVdrive';
    end
    if size(delta,1)~=1 %to do -- work out why modelVdotMax is sometimes sideways...
        delta=delta';
    end
    upperSEM=modelVdrive+delta/1.96;
    lowerSEM=modelVdrive-delta/1.96;
    if size(upperSEM,2)==1 %to do -- work out why modelVdotMax is sometimes sideways...
        upperSEM=upperSEM';
        lowerSEM=lowerSEM';
    end
    
    CI_parameters = nlparci(BestParametersDelay,RESIDUAL,'jacobian',JACOBIAN);
    CI_parameters = CI_parameters(:,2) - BestParametersDelay;
    LG1_CI = max(CI_parameters(1,:));
    tau_CI = max(CI_parameters(2,:));
    VRA_CI = max(CI_parameters(4,:));
    Delay_CI = max(CI_parameters(5,:));
end



%**************************************************************************
% Following the optimisation procedure the optimal values of LG_0, tau, Vdr[1],
% gamma and delay are now known. Equation S4 can now be used to calculate
% the dynamic loop gain parameters LG1, LG2, LGn and natural cycling
% frequency:
% LG(f)=-LG_0*e^(delay*i*2*pi*f)/(1+i*2*pi*f*tau)
%**************************************************************************


LG1=BestParameters(1); % The steady state loop gain
tau=BestParameters(2); % The estimated time-constant, Tau
gamma=BestParameters(4)/Veupnea; % The ventilatory response to arousal (VRA)
delay=BestParameters(5)*Ttot_median; % Convert the optimal delay from number of breaths back to a duration in seconds my multiplying by the patients median T_tot

tau2 = settings.tau2ontau1*tau;

Ff=[0.002:0.00001:1.09999 1.1.^(1:100)];
phase=pi-atan(2*pi*tau*Ff)-atan(2*pi*tau2*Ff)-2*pi*delay*Ff;
fn=interp1(phase,Ff,0,'spline'); % The natural cycling frequency
Tn=1/fn;
%fn = 1/Tn;
% tau=tan(pi-2*pi/Tn*delay)/(2*pi/Tn);
LG0 = LG1 * (1+(2*pi*tau/60)^2)^0.5 * (1+(2*pi*tau2/60)^2)^0.5;

% s = jw = j*2*pi*f
% LG = LG0 * exp(-s*delay) / (1+s*tau) / (1+s*tau2) ; % The loop gain at the natural frequency

LGn=abs(LG0*exp(-1i*2*pi*fn*delay)/(1+1i*2*pi*fn*tau)/(1+1i*2*pi*fn*tau2)/(1)); % The loop gain at the natural frequency
%LG1=abs(LG0*exp(-1i*2*pi*(1/60)*delay)/(1+1i*2*pi*(1/60)*tau)/(1)); % The loop gain at 1 cycles per minute
LG2=abs(LG0*exp(-1i*2*pi*(1/30)*delay)/(1+1i*2*pi*(1/30)*tau)/(1+1i*2*pi*(1/30)*tau2)/(1)); % The loop gain at 2 cycles per minute
%**************************************************************************
%**************************************************************************

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
LoopGainParameters=[LG0 tau delay LGn Tn LG1 LG2 gamma];
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%Vchem=Vchem+Veupnea;

end

