function tbl = FixMissingDataInTbl(tbl, cat, Gtest_train)

if length(unique(Gtest_train))<3 % we're missing a category in training data
    %keyboard
end

if length(unique(cat))<3 % we're missing a category in test data
    missingcats = setxor(unique(cat), [1 2 3]);
    if all(missingcats == [1,2])
        % move tbl data to col 3
        tbl2 = zeros(3);
        tbl2(:,3)=tbl;
        tbl = tbl2;
    elseif all(missingcats == [2,3])
        % set tbl data as col 1
        tbl2 = zeros(3);
        tbl2(:,1)=tbl;
        tbl = tbl2;
    elseif all(missingcats == [1,3])
        % set tbl data as col 2
        tbl2 = zeros(3);
        tbl2(:,2)=tbl;
        tbl = tbl2;
    elseif missingcats == 1
        % set tbl data as col 2 and 3
        tbl2 = zeros(3);
        tbl2(:,[2, 3])=tbl;
        tbl = tbl2;
    elseif missingcats == 2
        % set tbl data as col 1 and 3
        tbl2 = zeros(3);
        tbl2(:,[1, 3])=tbl;
        tbl = tbl2;
    elseif missingcats == 3
        % set tbl data as col 1 and 2
        tbl2 = zeros(3);
        tbl2(:,[1, 2])=tbl;
        tbl = tbl2;
    else
        keyboard
    end
end

% % first column are labels for the rows of tbl, which come from X (gtest)
% % second column are labels for the columns of tbl, which come from Y (cat)
% % get lbls from crosstab
% pp=cellfun(@isempty,lbls); %[i,j] = find(pp);
% find(pp(:,2)==1)