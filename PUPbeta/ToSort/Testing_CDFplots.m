%% 
% this is just some testing of CDF plots of Ttrans in VEVdrive data
%

InputDataCol = 74; % 74 is Ttrane_Ttot_ ^0.5 
   % could try other inputs...

%%
if 0 % group all pt data
    DataToPlot = Amatrix2(:, InputDataCol);
    [f1,x1] = ecdf(DataToPlot(Gtest_All>0.9).^2);
    [f2,x2] = ecdf(DataToPlot((Gtest_All>0.5)&(Gtest_All<0.9)).^2);
    [f3,x3] = ecdf(DataToPlot(Gtest_All<0.5).^2);
    figure(1); clf(figure(1));
    plot(x1, f1, 'g-'); hold on;
    plot(x2, f2, 'b-');
    plot(x3, f3, 'r-');
end

%% do for indiv pts
p1 =[]; p2 =[]; p3 = [];

for subj=1:54 % set pt num or all, no data for PT=1
    if ~ismember(subj, PT_list) %PT_list)
        continue
    end
    %subj
    Isubj=(PtData.PT==subj);
    Gtest_pt = Gtest_All(Isubj);    
    DataToPlot = Amatrix2(Isubj, InputDataCol);
    
    [f1,x1] = ecdf(DataToPlot(Gtest_pt>0.9).^2); 
    [f2,x2] = ecdf(DataToPlot((Gtest_pt>0.5)&(Gtest_pt<0.9)).^2);
    [f3,x3] = ecdf(DataToPlot(Gtest_pt<0.5).^2);
    p1 = [p1; prctile(x1,[5 10 25 50 75])];
    p2 = [p2; prctile(x2,[5 10 25 50 75])];
    p3 = [p3; prctile(x3,[5 10 25 50 75])];
    
    figure(2); clf(figure(2)); 
    plot(x1, f1, 'g-'); hold on;
    plot(x2, f2, 'b-');
    plot(x3, f3, 'r-');

end

%%
if 1 % these pts
    %pts = [26 28 2 4 17 38 18 29 30 34 35 47]; % example patients
    pts = [26 28 2 4 17 38 18 29 30 34 35 47 3 5 6 8 10 11 14 19 22 23 24 27 39 40 41 43 46 50 53 54]; %
    %pts = [9 16 20 21 33 45 44]; % these guys actually went opposite to expected
    pts_ = [];
    for n=1:length(pts)
        pt = find(PT_list==pts(n));
        pts_ = [pts_, pt];
    end
else
    pts_ = 1:41; % all
end

tile = 2; % 1=5%, 2=10%, 3=25%, 4=50%, 5=75%
figure(3); clf(figure(3));
boxplot([p1(pts_,tile),p2(pts_,tile),p3(pts_,tile)]);
