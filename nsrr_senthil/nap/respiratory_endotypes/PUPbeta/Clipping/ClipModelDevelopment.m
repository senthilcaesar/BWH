addpath(genpath(('C:\Users\SOpdeBeeck\Dropbox\PUPbeta_git\PUPbeta20190629')));

dir = 'D:\Clipping algorithm\Analyses\Analyzed\';
savename = 'SaraODB_OAT';

M=36;
for m=1:M
    disp(num2str(m));
    load([dir savename '_' num2str(m)],'BreathDataTable');
    [~,BreathDataTableE_]=GetNonOvlappedVE(BreathDataTable);
    BreathDataTableE_.Subj = m*ones(size(BreathDataTableE_,1),1);
    if ~exist('DataTblAll')
        DataTblAll = BreathDataTableE_;
    else
        DataTblAll = [DataTblAll;BreathDataTableE_]; 
    end
end
%%
DataTblAllOrig = DataTblAll;
%%
DataTblAll = DataTblAll(:,{'FclippedUpper','FclippedLower','VTiCorrection','VTeCorrection','Subj','Time_start','Time_mid','Time_end'});
%%
save DataTblAll DataTblAll
%%
figure(998); clf(998); set(gcf,'color',[1 1 1]);
           scatter(DataTblAll.FclippedUpper,DataTblAll.VTiCorrection,4,[0.8 0.2 0.1],'filled','markerfacealpha',0.3)
           
           axis([0 1 0.9 2]);
           hold on
           scatter(DataTblAll.FclippedLower,DataTblAll.VTeCorrection,4,[0.2 0.8 0.1],'filled','markerfacealpha',0.3)
           
           %plot(DataTblAll.FclippedLower,DataTblAll.VTeCorrection,'.')
           
           Fline = 0:0.01:0.9;
           
           Yline = 1./( cos(0.5*pi()*Fline).^0.33);
          
           
           f_exp2=fit(DataTblAll.FclippedLower(isnan(DataTblAll.VTeCorrection)==0),...
               DataTblAll.VTeCorrection(isnan(DataTblAll.VTeCorrection)==0),'poly2');
           Yline_Exp = f_exp2.p1*(Fline).^2+f_exp2.p2*(Fline)+f_exp2.p3;
           
           
           f_insp2=fit(DataTblAll.FclippedUpper(isnan(DataTblAll.VTiCorrection)==0),...
               DataTblAll.VTiCorrection(isnan(DataTblAll.VTiCorrection)==0),'poly2');
           Yline_Insp = f_insp2.p1*(Fline).^2+f_insp2.p2*(Fline)+f_insp2.p3;
           
           plot(Fline,Yline_Exp,'b')
           hold on
           plot(Fline,Yline_Insp,'m')
  
           hold off
           box off
%% Bins
Edges = [0 0.01 0.1:0.1:0.9];
EdgesN = length(Edges);
IDa = [sum(DataTblAll.FclippedUpper>Edges,2)];
for i=1:EdgesN
    I = IDa==i;
    Ybin(i) = nanmedian(DataTblAll.VTiCorrection(I));
    Xbin(i) = nanmedian(DataTblAll.FclippedUpper(I));
    Nbin(i) = sum(I);
end

%% Bins
Edges = [0 0.01 0.1:0.1:0.9];
EdgesN = length(Edges);
IDaE = [sum(DataTblAll.FclippedLower>Edges,2)];
for i=1:EdgesN
    I = IDaE==i;
    YbinE(i) = nanmedian(DataTblAll.VTeCorrection(I));
    XbinE(i) = nanmedian(DataTblAll.FclippedLower(I));
    NbinE(i) = sum(I);
end

%% Figure inspiration
figure(999); clf(999); set(gcf,'color',[1 1 1]);
I = ~isnan(DataTblAll.VTiCorrection) & DataTblAll.FclippedUpper>0.01 & DataTblAll.VTiCorrection>0.5 & DataTblAll.VTiCorrection<10;

Ftrans = @(x) (x-0.5).^0.1 -0.5^0.1;
FtransInv = @(x) ((x+0.5^0.1).^10+0.5);

scatter((DataTblAll.FclippedUpper(I)),Ftrans(DataTblAll.VTiCorrection(I)),6,[0.1 0.1 0.8],'filled','markerfacealpha',0.1)

Fline = [0:0.01:0.95];
hold on
plot(Xbin,Ftrans(Ybin),'g.-')

if 0
Yline = (1./(cos(0.5*pi()*Fline).^0.33)-1)*0.13;
plot(Fline,Yline,'r--')          
end

fiteq = 'a*tan(x*1.5)'; start = [0.02];
[f_insp2,gof]=fit(DataTblAll.FclippedUpper(I),Ftrans(DataTblAll.VTiCorrection(I)),fiteq,'start',start);
Yline_Insp2 = feval(f_insp2,Fline);


hold on
plot(Fline,Yline_Insp2,'r','linewidth',2)
figure(999);

if 1
set(gca,'ytick',Ftrans([0.5 0.67 0.75 1 1.5 2 3 4]))
set(gca,'yticklabels',FtransInv(get(gca,'ytick')));
end
ylim([Ftrans([0.8 4])])

%
figure(998); clf(998); set(gcf,'color',[1 1 1]);
scatter((DataTblAll.FclippedUpper(I)),(DataTblAll.VTiCorrection(I)),6,[0.1 0.1 0.8],'filled','markerfacealpha',0.1)
hold on
plot(Fline,FtransInv(Yline_Insp2),'r','linewidth',2)

VTiCorrectionF = @(x) (((0.02063*tan(x*1.5))+0.5^0.1).^10+0.5);
%plot(Fline,VTiCorrectionF(Fline),'k--','linewidth',2);


%% Figure exspiration
figure(997); clf(997); set(gcf,'color',[1 1 1]);
I = ~isnan(DataTblAll.VTeCorrection) & DataTblAll.FclippedLower>0.01 & DataTblAll.VTeCorrection>0.5 & DataTblAll.VTeCorrection<10;

Ftrans = @(x) (x-0.5).^0.1 -0.5^0.1;
FtransInv = @(x) ((x+0.5^0.1).^10+0.5);

scatter((DataTblAll.FclippedLower(I)),Ftrans(DataTblAll.VTeCorrection(I)),6,[0.1 0.1 0.8],'filled','markerfacealpha',0.1)

Fline = [0:0.01:0.95];
hold on
plot(XbinE,Ftrans(YbinE),'g.-')

if 0
Yline = (1./(cos(0.5*pi()*Fline).^0.33)-1)*0.13;
plot(Fline,Yline,'r--')          
end

fiteq = 'a*tan(x*1.2)'; start = [0.02];
[f_insp2,gof]=fit(DataTblAll.FclippedLower(I),Ftrans(DataTblAll.VTeCorrection(I)),fiteq,'start',start);
Yline_Insp2 = feval(f_insp2,Fline);


hold on
plot(Fline,Yline_Insp2,'r','linewidth',2)

if 1
set(gca,'ytick',Ftrans([0.5 0.67 0.75 1 1.5 2 3 4]))
set(gca,'yticklabels',FtransInv(get(gca,'ytick')));
end
ylim([Ftrans([0.8 4])])

%
figure(996); clf(996); set(gcf,'color',[1 1 1]);
scatter((DataTblAll.FclippedLower(I)),(DataTblAll.VTeCorrection(I)),6,[0.1 0.1 0.8],'filled','markerfacealpha',0.1)
hold on
plot(Fline,FtransInv(Yline_Insp2),'r','linewidth',2);

VTeCorrectionF = @(x) (((0.0527*tan(x*1.2))+0.5^0.1).^10+0.5);
%plot(Fline,VTeCorrectionF(Fline),'k--','linewidth',2);

%%

save ClippingCorrectionF VTiCorrectionF VTeCorrectionF