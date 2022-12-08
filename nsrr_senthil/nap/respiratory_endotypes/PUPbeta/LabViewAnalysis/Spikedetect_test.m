function [respwav]=Spikedetect_test(respwav,fs)

Irmv=zeros(length(respwav),1);
n1=fs; n2=fs/2;
xmax=movmax(respwav,[n1,n1]);
xmax=movmin(xmax,[n2,n2]);
% figure; plot(timewav,respwav); hold on;
% plot(timewav,xmax,'r');

n1=fs; n2=fs;
xmax2=movmax(respwav,[n1,n1]);
xmax2=movmin(xmax2,[n2,n2]);
% figure; plot(timewav,respwav); hold on;
% plot(timewav,xmax2,'r');
%
xmax_tmp=(xmax-xmax2);  figure(88); clf(figure(88)); plot(respwav); hold on; plot(xmax_tmp,'r');
sig=diff(xmax_tmp);
hold on; plot(sig,'g');

tmp_change = mode(sig);
tmp_xmax_tmp = mode(xmax_tmp);

index_xmap_tmp = find(xmax_tmp>1.5);

respwav(index_xmap_tmp(1:end)) = NaN;

%     I=[];flowtemp=[];
%     for j=1:length(MinIdx)
%         mxri=MaxIdx(j)+10;
%         mnri=MinIdx(j)+15;
%         if mnri>mxri, ri=mxri; else, ri=mnri; end
%         flowtemp=respwav(MinIdx(j):ri);
%         [mx1,Itemp]=max(flowtemp);
%         if mx1>prctile(flowtemp,75),I=[I;MinIdx(j)+Itemp-1];end
%         flowtemp=[];
%         if ri==mnri % larger duration b/w min and max index indicate double spikes
%             
%             flowtemp=respwav(ri+1:mxri);
%             [mx2,Itemp2]=max(flowtemp);
%             if mx2>prctile(flowtemp,75),I=[I;ri+1+Itemp2-1];end
%             flowtemp=[];
%             
%         end
%     end
%     
%     I(respwav(I)<(prctile(respwav(I),5)))=[]; % remove small spikes detected
%     Idiff=diff(I);
%     I(Idiff==1)=[]; % merge close spikes if any
%     
%     
%     MinIdxn=[]; MaxIdxn=[];
%     for i=1:length(I) % finding the start and end of diff signal frm which spikes were detected
%         fmn=find(MinIdx<=I(i),1,'last');
%         MinIdxn=[MinIdxn; MinIdx(fmn)];
%         fmx=find(MaxIdx>=MinIdx(fmn),1,'first');
%         MaxIdxn=[MaxIdxn; MaxIdx(fmx)];
%     end
%     hold on; plot(timewav(I),respwav(I),'*k')
%     %makes spikes wider by Fextrachop (upper limit maxextrachop
%     %sec)
%     I1=MinIdxn;
%     I2= MaxIdxn;
%     lengthI=(I2-I1)*(1/fs);
%     lengthIi=(I2-I1);
%     Fextrachop=1.5;
%     maxextrachop=0.1; %sec
%     extrachopi=floor(lengthIi*Fextrachop);
%     maxextrachopi=round(maxextrachop*fs);
%     extrachopi(extrachopi>maxextrachopi)=maxextrachopi;
%     
%     for i=1:length(I)
%         li=I1(i)+1-extrachopi(i);
%         if li<I(i)-10
%             li=(I(i)-10);
%         end
%         if li<1
%             li=1;
%         end
%         ri=I2(i)+extrachopi(i);
%         if ri>I(i)+10
%             ri=I(i)+10;
%         end
%         if ri>length(respwav)
%             ri=length(respwav);
%         end
%         Irmv(li:ri)=1;
%     end
%     
%     hold('on');
%     respwav2=respwav;
%     timewav2=timewav;
%     delval=respwav2(Irmv==1);
%     
%     % checking to see if there are any more spikes
%     
%     rprc=respwav2(respwav2>=prctile(delval,85));
%     tprc=timewav2(respwav2>=prctile(delval,85));
%     if ~isempty(rprc)
%         hold on; plot(tprc,rprc,'*r');
%         
%         rprc(round(diff(tprc))<1)=[];
%         tprc(round(diff(tprc))<1)=[];
%         [~,idxs] = intersect(timewav2,tprc,'stable');
%         
%         for i=1:length(idxs)
%             respwav2(idxs(i)-16:idxs(i)+16)=NaN;
%             timewav2(idxs(i)-16:idxs(i)+16)=NaN;
%             Irmv(idxs(i)-16:idxs(i)+16)=1;
%         end
%     end
%     hold on; plot(timewav2,respwav2);
%     
%     respwav2(Irmv==1)=NaN; % remove spikes
%     plot(timewav2,respwav2);
%     
%     timewavtemp=timewav(Irmv==0); % keeping only the non spiky data
%     respwavtemp=respwav(Irmv==0);
%     
%     timewavtemp2=timewavtemp(~isnan(respwavtemp));
%     respwavtemp2=respwavtemp(~isnan(respwavtemp));
%     respwav3=interp1(timewavtemp2,respwavtemp2,timewav,'pchip');
%     hold('on');plot(timewav,respwav3);
%     respwav=respwav3;
%     lrespwavtmp=length(respwavtemp2);
%     
%   
% else
%     lrespwavtmp=length(respwav);
%     respwav=respwav;
%   
%     
% end
end