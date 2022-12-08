function [SaO2XHz3]=SpO2ArtifactReject(SaO2XHz,dt)

SaO2XHz(isnan(SaO2XHz))=0; % remove NaN for filter 

filter_HFcutoff_butter1=1; %1 Hz smoothing filter
filter_order1 = 2;
[B_butter1,A_butter1] = butter(filter_order1,[filter_HFcutoff_butter1]/(1/dt/2),'low');
SaO2XHz=filter(B_butter1,A_butter1,SaO2XHz);

SaO2XHz = round(SaO2XHz);

% Auto SpO2 artefact removal
SpO2artthres1=40; %set this
proximity_thres=10; %set this
lower_thres_prctile=0; %set this

SaO2XHz1=SaO2XHz;
for i=1:length(SaO2XHz)
    if SaO2XHz(i)<SpO2artthres1
        SaO2XHz1(i)=NaN;
    end
end

SaO2XHz2=SaO2XHz;
SpO2neighborhoodD_dt=2; %set this
SpO2neighborhoodD_i=round(SpO2neighborhoodD_dt/dt);
for i=(SpO2neighborhoodD_i+1):(length(SaO2XHz1)-SpO2neighborhoodD_i)
    if isnan(sum(SaO2XHz1(i-SpO2neighborhoodD_i:i+SpO2neighborhoodD_i)))
        SaO2XHz2(i)=NaN;
    end
end

SaO2XHz2(SaO2XHz2>100&SaO2XHz2<101)=100;
SaO2XHz2(SaO2XHz2>101)=NaN;
I=~isnan(SaO2XHz2);
SaO2info1=zeros(1,length(SaO2XHz2))+1000000;
temp=0;
startcounting=0;
for i=2:length(SaO2info1)
    if (I(i-1)==0&&I(i)==1)&&~startcounting
        temp=0;
        SaO2info1(i)=temp;
        startcounting=1;
    elseif (I(i-1)==1&&I(i)==1)&&startcounting
        temp=temp+dt;
        SaO2info1(i)=temp;
    elseif I(i)==0
        startcounting=0;
        temp=0;
        SaO2info1(i)=temp;
    end
end
SaO2info2=zeros(1,length(SaO2XHz2))+1000000;
temp=0;
startcounting=0;
for i=(length(SaO2info2)-1):-1:1
    if (I(i+1)==0&&I(i)==1)&&~startcounting
        temp=0;
        SaO2info2(i)=temp;
        startcounting=1;
    elseif (I(i+1)==1&&I(i)==1)&&startcounting
        temp=temp+dt;
        SaO2info2(i)=temp;
    elseif I(i)==0
        startcounting=0;
        temp=0;
        SaO2info2(i)=temp;
    end
end

SaO2XHz3=SaO2XHz2;
I1=SaO2info1(:)<proximity_thres|SaO2info2(:)<proximity_thres;
lower_thres=prctile(SaO2XHz2(~I1),lower_thres_prctile);
I2=SaO2XHz2<lower_thres;
I3=I2(:)&I1(:);
SaO2XHz3(I3)=NaN;