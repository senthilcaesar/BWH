function [selectsnoreidx,batch] = SelectSnores(EnvelpData,Sites,BreathDataTable,Freq,analyzedIdx,SnoreCounts)

selectsnoreidx = false(size(EnvelpData,1),1);
SubjectIDBr = BreathDataTable.SubjectID(analyzedIdx);
SubUniq = unique(SubjectIDBr);
SiteUniq = unique(Sites);
% FeatureTableSub = nan(length(SubUniq)*length(SiteUniq),size(FeatureTableBr,2)-1);
SitesSub = repmat(SiteUniq,[length(SubUniq),1]);
% selectsnoreidx = [];
batch = nan(size(EnvelpData,1),1);
batchnum = 0;

for ii = 1:length(SubUniq)
    subidx = strcmp(SubjectIDBr,SubUniq{ii});
    disp([SubUniq{ii},': ',num2str(sum(subidx)), ' snores'])
    for jj = 1:length(SiteUniq)
        siteidx = Sites==SiteUniq(jj);
        subsiteidx = subidx & siteidx;
        disp([SubUniq{ii},', ',char(SiteUniq(jj)),': ',num2str(sum(subsiteidx)), ' snores'])
        
        numsnoreskeep = SnoreCounts.SnoresPerSub(strcmp(char(SiteUniq(jj)),SnoreCounts.Properties.RowNames));
        
        if sum(subsiteidx) < 1
            continue
        elseif sum(subsiteidx) < numsnoreskeep
            numsnorespick = sum(subsiteidx);
            batchnum = batchnum + 1;
        elseif sum(subsiteidx) >= numsnoreskeep
            numsnorespick = numsnoreskeep;
            batchnum = batchnum + 1;
        end
        
        EnvelpData_ = EnvelpData(analyzedIdx,:);
        EnvelpSiteSub = EnvelpData_(subsiteidx,:);
        
        medianEnvelp = nanmedian(EnvelpSiteSub,1);
        envRMSEall = sqrt(nansum((EnvelpData - medianEnvelp).^2,2)./size(EnvelpData,2));
        envRMSEsub = sqrt(nansum((EnvelpSiteSub - medianEnvelp).^2,2)./size(EnvelpData,2));
        
        minEnvRMSE = sort(envRMSEsub,'ascend');
        minEnvRMSE_10 = minEnvRMSE(1:numsnorespick);
        selectidxtemp = ismember(envRMSEall, minEnvRMSE_10);
        selectsnoreidx = selectsnoreidx + selectidxtemp;
        batch(selectidxtemp) = batchnum;
%         selectsnoreidx = [selectsnoreidx; selectidxtemp];

        % plot
        if 0
            figure(10); hold off
            plot(Freq, medianEnvelp,'k','LineWidth', 2); hold on
            selectidxnum = find(selectidxtemp);
            subsiteidxnum = find(subsiteidx);
            for kk = 1:sum(subsiteidx)
                if sum(ismember(subsiteidxnum(kk),selectidxnum))==1
                    fig1 = plot(Freq, EnvelpData(subsiteidxnum(kk),:),'k');
                    fig1.Color = [fig1.Color 0.5];
                    xlim([0 3.5])
                    ylim([0 40])
                else
                    fig2 = plot(Freq, EnvelpData(subsiteidxnum(kk),:),'r');
                    fig2.Color = [fig2.Color 0.1];
                    xlim([0 3.5])
                    ylim([0 40])
                end
            end
        end
        
    end
end