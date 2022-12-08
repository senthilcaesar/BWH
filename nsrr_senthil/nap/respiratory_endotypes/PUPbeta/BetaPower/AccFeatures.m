function TAccFeatures = AccFeatures(WSpredlogit,Exclude)
global settings

%% AccFeatures
%mdlAcc
%clear maxcorr r2 h1 bimodalitycoef mean_ SD_ diffhistmin

%WSpredlogit(Iw,i), Exclude(:,i), j is Neegs
plotson=0;  % DLM changed 1 to 0, on 20200331 1536 GMT+10
j = size(WSpredlogit,2);

        for i=1:j
            if plotson
                subplot(1,j,i)
            end
            
            dStep=0.25;
            Centers=-10:dStep:10;
            Edges=(Centers(1)-dStep/2):dStep:(Centers(end)+dStep/2);
            
            [h1(i,:),edges] = histcounts(WSpredlogit(~Exclude(:,i),i),Edges);
            h1(i,:) = h1(i,:)/sum(h1(i,:));
            
            if plotson
                bar(Centers,h1(i,:),'EdgeAlpha',0,'FaceAlpha',0.3,'BarWidth',1);
                hold('on');
            end
            
            Iw = WSpredlogit(:,i)>0 & ~Exclude(:,i);
            Is = WSpredlogit(:,i)<=0 & ~Exclude(:,i);

            %%%
            mean_(i,1) = nanmean(WSpredlogit(~Exclude(:,i),i));
            SD_(i,1) = nanstd(WSpredlogit(~Exclude(:,i),i));
            
            try
            ft = fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2)');
            ftoptions = fitoptions(ft);
            ftoptions.StartPoint = [0.1 0.1 5 -5 3 3];
            ftoptions.Lower = [0 0 0  -10  0 0];
            ftoptions.Upper = [1 1 10  0 10 10];
            [fitobject,gof] = fit(Centers',h1(i,:)',ft,ftoptions);
            r2(i,1) = gof.rsquare;
            if plotson
                plot(fitobject)
            end
            catch me
               r2(i,1) = NaN;
            end
            
            %%%
            bimodalitycoef(i,1) = (skewness(WSpredlogit(~Exclude(:,i),i))^2 + 1)/kurtosis(WSpredlogit(~Exclude(:,i),i));
            
            %%%
            noncentralness(i,1) = (nanmean((WSpredlogit(~Exclude(:,i),i)).^2))^0.5;
            
            %%%
            temp = WSpredlogit(~Exclude(:,i),i);
            temp(temp<-10)=-10; temp(temp>10)=10;
            temp = min([abs(temp-10) abs(temp+10)]')';
            nonperipheralness(i,1) = (nanmean(temp.^2))^0.5;
            
%             if plotson
%                 title(num2str(acc(i,1),2),'FontWeight','normal');
%             end
            
            
            %%% correlations with others
            if j==1 %only one signal, thus can not compare to others
                maxcorr(i,1)=0;
            else
            clear corrtemp
            temp = WSpredlogit(:,i);
            temp(Exclude(:,i))=NaN;
            for x=1:j
                corrtemp(x)=NaN; %default
                if x==i
                    continue                    
                end
                
                temp2 = WSpredlogit(:,x);
                temp2(Exclude(:,x))=NaN;
                isnans = isnan(temp)|isnan(temp2);
                if sum(isnans)/length(isnans)==1
                    continue
                else
                    try
                    corrtemp(x) = corr(temp(~isnans),temp2(~isnans));
                    catch me
                    me.message
                    end
                end
                %error message: "Requires a data matrix X" - when?
                
            end
            tempmaxcorrtemp = max(corrtemp);
            if isnan(tempmaxcorrtemp)
                tempmaxcorrtemp=0;
            end
            maxcorr(i,1)=tempmaxcorrtemp; %%%
            end
        end
        
        
        for i=1:j %%%
            range=1:j; 
             if isfield(settings,'AllowOneEEG') && settings.AllowOneEEG == 1 % DLM testing for limited EEG channel data
                 % don't remove range(i)                
             else               
                range(i)=[]; % normal operation  
             end
            %diffhistmin(i,1) = sum(min(abs(h1(range,:) - h1(i,:))));
            diffhistmin(i,1) = min(sum(abs(h1(range,:) - h1(i,:)),2));
        end
        
        TAccFeatures = table(maxcorr,bimodalitycoef,r2,diffhistmin,mean_,SD_,noncentralness,nonperipheralness);
        