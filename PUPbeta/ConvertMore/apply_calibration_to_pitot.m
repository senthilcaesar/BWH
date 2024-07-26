function CalibratedFlow = apply_calibration_to_pitot(Flow, fname)

        % Here, we are attempting to apply the Spike formula 
        % e.g. "Poly(Hwr(Ch(PTch)-1*0.491059),0.0,2.874796,9.098049,-11.245381)+Poly((Hwr((Ch(PTch)-1*0.491059)*-1))*-1,0.0,4.241441,-5.333880,-8.465138)"
        % to the pitot flow signal, in order to convert the imported signal to calibrated flow
        % CAUTION: if successful, this will overwrite the exisitng flow signal
        
%% Calibrated flow when Pitot signal is adjusted by unique formula (provided in spreadsheet per ID)
% Basic process is:
% 1. get the uncalibrated pitot signal (passed in as Flow)
% 2. read in corresponding formula (read Cal spreadsheet, find current ID)
% 3. apply calibration formula to pitot signal

% read the formulas spreadsheet
Formulas = readtable('Pitot_Cal_Equations.xlsx');
% this is a spreadsheet with rows per study ID (matching that in the AMasterSpreadsheet)
% column A is Study ID (currently with Header row, "Study_ID")
% column B is the Formula (header row contains, "Pitot_Cal_Equations")

% find the appropriate formula for this fname
ID = fname(1:end-4);
indx = find(strcmp(Formulas.Study_ID, ID));
CalEq = Formulas.Pitot_Cal_Equations(indx);
% {'Poly(Hwr(Ch(PTch)-1*0.491059),0.0,2.874796,9.098049,-11.245381)+Poly((Hwr((Ch(PTch)-1*0.491059)*-1))*-1,0.0,4.241441,-5.333880,-8.465138)'}

%% Spike
% Poly(x,L)  Replace x with a polynomial in x of order 1 to 5; L is a list of 1 to 6 coefficients. 
%            Use this to apply a non-linear calibration to a signal. 
%            For example: Poly(x,a,b,c,d) = a + b*x + c*x*x + d*x*x*x
%
% Hwr(x)     Half wave rectify x. Negative values are replaced by zeros. This is faster than Max(x, 0).

%% matlab
% y = polyval(p,x) evaluates the polynomial p at each point in x. 
% The argument p is a vector of length n+1 whose elements are the coefficients (in descending powers) of an nth-degree polynomial:
% p(x)=p1.x.x + p2.x + ... + pn.x + pn+1.

%% split the formula into the two main parts
EQparts = strsplit(string(CalEq), '+');

%% process first part (second part is same method)
EQparts(1) = strrep(EQparts(1),')','');
%strip first N=16 chars  [ Poly(Hwr(Ch(PTch) ],
EQ1 = char(EQparts(1)); EQ1(1:16)=[]; EQ1Terms = strsplit(EQ1, ','); 
% split into offset and individual terms
offset = eval(cell2mat(EQ1Terms(1))); % signal offset
a = str2double(cell2mat(EQ1Terms(2))); % a
b = str2double(cell2mat(EQ1Terms(3))); % b
c = str2double(cell2mat(EQ1Terms(4))); % c
d = str2double(cell2mat(EQ1Terms(5))); % d

FlowAdj = Flow+offset;                      % apply signal offset
FlowPos = max(FlowAdj, 0);                  % get the positive component of the pitot signal (HWR in Spike, negative values are set to 0)
FlowPosCal = polyval([d,c,b,a], FlowPos);   % apply part 1 polynomial

%% second part - not strictly fully automated... 
% in spike, it's flipping the signal to use hwr while getting the negative portion, then reflipping to push it back negative
EQparts(2) = strrep(EQparts(2),')','');
%strip first N=16 chars  [ Poly(Hwr(Ch(PTch) ],
EQ2 = char(EQparts(2)); EQ2(1:18)=[]; EQ2Terms = strsplit(EQ2, ','); 
% split into offset and individual terms
offset = eval(cell2mat(EQ2Terms(1))); % signal offset
a = str2double(cell2mat(EQ2Terms(2))); % a
b = str2double(cell2mat(EQ2Terms(3))); % b
c = str2double(cell2mat(EQ2Terms(4))); % c
d = str2double(cell2mat(EQ2Terms(5))); % d

FlowAdj = Flow+offset;                      % apply signal offset (in test case this offset was the same as part one, but applying here for the possibility that they could be different)
FlowNeg = min(FlowAdj, 0);                  % get the negative component of the pitot signal
FlowNegCal = polyval([d,c,b,a], FlowNeg);   % apply part 2 polynomial

%% re-assemble the positive and negative components
CalibratedFlow = FlowPosCal + FlowNegCal;








