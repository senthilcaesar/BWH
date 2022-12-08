function [kappahistory_pt,B4,currentlist] = SequentialForwardSelection(Amatrix2_train,Gtest_train,Ftrs)
global settings
addthreshold = 0.005;
kappahistory_pt = []; % reset kappa history within pt

%% sequential forward floating search
% Ftrs = Ftr_indx_flow;
%if settings.findfirstftr
    % do one loop through (all this is just to get the starting ftr)
    % side note, this may be helpful when finding 'super-features'
    progressbar('','FirstFtr','');
    progressbar([],2/10,[]);
    kappa_array = nan(length(Ftrs),1);
    for fwdcycle = 1:length(kappa_array)
        progressbar([],[],fwdcycle/length(kappa_array));
        templist = Ftrs(fwdcycle);
        if settings.verbose; str=['Testing ftr #', num2str(templist)]; disp(str); end % include indicator of loop
        [B1,dev,stats]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'nominal');
        [pihat] = mnrval(B1,Amatrix2_train(:,templist), 'Model', 'nominal');
        [prob,cat] = max(pihat(:,:),[],2);
        tbl = confusionmat(grp2idx(Gtest_train),cat);
%         try % if it's not 3x3, then we've got a problem
%             tbl=tbl_.*ones(3,4);
%         catch % this happens when not all cats are in pred or train,
%             % so, work out what we've got, and what we are missing,
%             % and fill it in so that kappa can work at least a bit
%             tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
%         end
        [kappa_array(fwdcycle), ~] = kappaDLMmod(tbl);
    end
%    save('FirstFtrKappa','kappa_array');
%else
%    load('FirstFtrKappa.mat');
%end % find first feature loop

[currentkappa, newindx] = max(kappa_array);
currentlist = Ftrs(newindx);
if settings.verbose; str=['Starting ftr #', num2str(currentlist)]; disp(str); end
% add this kappa to the record
kappahistory_pt = [kappahistory_pt, currentkappa];

try
    progress = 1;
    while progress % while we are making progress
        progressbar([],(length(currentlist)+2)/10,[]);
        % try adding a single new ftr
        %   repeat for N of remaining list
        %  set up the list of available ftrs, (i.e. all those not in the current list)
        availftrs = Ftrs; indxout = [];
        for ind = 1:length(currentlist) % there must be a better way...
            loc = find(Ftrs==currentlist(ind));
            indxout = [indxout, loc];
        end
        availftrs(indxout)=[];

        kappa_array = nan(length(availftrs),1);
        progressbar('',['Test adding ftrs, loop ', num2str(length(currentlist))],'');
        for fwdcycle = 1:length(kappa_array)
            progressbar([],[],fwdcycle/length(kappa_array));

            % cycle through current list plus one from remaining list
            templist = [currentlist, availftrs(fwdcycle)];
            if settings.verbose; str=['Testing ftrs + ', num2str(templist)]; disp(str); end

            [B2,dev,stats]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'nominal');
            [pihat] = mnrval(B2,Amatrix2_train(:,templist), 'Model', 'nominal');
            [prob,cat] = max(pihat(:,:),[],2);
            tbl = confusionmat(grp2idx(Gtest_train),cat);
%             try % if it's not 3x3, then we've got a problem
%                 tbl=tbl_.*ones(3);
%             catch % this happens when not all cats are in pred or train, this happens in pt 12,
%                 % so, work out what we've got, and what we are missing,
%                 % and fill it in so that kappa can work at least a bit
%                 tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
%             end
            [kappa_array(fwdcycle), ~] = kappaDLMmod(tbl); % get cohens kappa valu
        end
        % find best new combination
        [newkappa, newindx] = max(kappa_array);
        newlist = [currentlist, availftrs(newindx)];
        % if new kappa better than current kappa, and
        %   add improvement delta > add threshold, then
        %  set as current, else add progress = 0
        adddelta = newkappa - currentkappa;
        if newkappa > currentkappa && adddelta > addthreshold
            currentkappa = newkappa;
            currentlist = newlist;
            addprogress = 1;
            if settings.verbose; str=['Adding ftr ', num2str(availftrs(newindx))]; disp(str); end
        else
            addprogress = 0;
        end

        % try removing a single feature to see if this improves things
        %  repeat for N of current list
        % note, there is actually no point doing this if currentlist is only two long,
        %  because we've already just tried all the individual features
        if length(currentlist) > 4
            kappa_array = nan(length(currentlist),1);
            progressbar('',['Test removing ftrs, loop ', num2str(length(currentlist))],'');
            for revcycle = 1:length(kappa_array)
                progressbar([],[],revcycle/length(kappa_array));

                % cycle through current list missing one feature
                templist = currentlist; templist(revcycle) = [];
                if settings.verbose; str=['Testing ftrs - ', num2str(templist)]; disp(str); end

                [B3]=mnrfit(Amatrix2_train(:,templist),Gtest_train, 'Model', 'nominal');
                [pihat] = mnrval(B3,Amatrix2_train(:,templist), 'Model', 'nominal');
                [prob,cat] = max(pihat(:,:),[],2);
                tbl = confusionmat(grp2idx(Gtest_train),cat);
%                 try % if it's not 3x3, then we've got a problem
%                     tbl=tbl_.*ones(3);
%                 catch % this happens when not all cats are in pred or train, this happens in pt 12,
%                     % so, work out what we've got, and what we are missing,
%                     % and fill it in so that kappa can work at least a bit
%                     tbl=FixMissingDataInTbl(tbl_,cat,Gtest_train);
%                 end
                [kappa_array(revcycle), ~] = kappaDLMmod(tbl); % get cohens kappa value
            end

            % find best new combination
            [newkappa, newindx] = max(kappa_array);
            newlist = currentlist; newlist(newindx)=[];
            % if new kappa better than current kappa, and
            %   remove improvement delta > remove threshold, then
            %  set as current, else remove progress = 0
            removedelta = newkappa - currentkappa;
            if newkappa > currentkappa && removedelta > removethreshold
                currentkappa = newkappa;
                currentlist = newlist;
                removeprogress = 1;
                if settings.verbose; str=['Removing ftr ', num2str(newlist(newindx))]; disp(str); end
            else
                removeprogress = 0;
            end
        else
            removeprogress = 0; %
        end

        % if add progress = 0 and remove progress = 0
        if addprogress == 0 && removeprogress == 0
            progress = 0; % then progress = 0, end process.
        else
            kappahistory_pt = [kappahistory_pt, currentkappa];
        end
    end % while progress within this subj
catch SFFS_error
    disp(SFFS_error.getReport)
end
progressbar('','PtModel','');
progressbar([],9/10,[]);

%% make one final model for this pt
[B4]=mnrfit(Amatrix2_train(:,currentlist),Gtest_train, 'Model', 'nominal');

end