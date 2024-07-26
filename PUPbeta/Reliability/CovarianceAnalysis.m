
Trait2P = SummaryAnalysisTable_2{crit==1,6}; % trait meeting criteria--even windows
Trait1P = SummaryAnalysisTable_1{crit==1,6}; % trait meeting criteria---odd window

Trait2A = SummaryAnalysisTable_2{crit==1,7}; % trait meeting criteria--even windows
Trait1A = SummaryAnalysisTable_1{crit==1,7}; % trait meeting criteria---odd window

Trait2C = SummaryAnalysisTable_2{crit==1,8}; % trait meeting criteria--even windows
Trait1C = SummaryAnalysisTable_1{crit==1,8}; % trait meeting criteria---odd window

%temp = Trait2C == Trait2A - Trait2P;

I = ~isnan(Trait2P) &  ~isnan(Trait2A);
C = cov(Trait2P(I),Trait2A(I))
C(2,1)


A = 2*randn(50,1);
B = 1*randn(50,1);
B = B+0.0*A;
corr(A,B)
corr(Trait2P(I),Trait2A(I));

C = cov(A,B)
[nanvar(A) nanvar(B) nanvar(A)+nanvar(B)-2*C(2,1) nanvar(A-B)]

[nanvar(Trait2A(I)) nanvar(Trait2A(I)) nanvar(Trait2P(I))+nanvar(Trait2A(I))-2*C(2,1) nanvar(Trait2C(I))]

mdlErrKeep2{1}



