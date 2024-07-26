addpath('/Users/sq566/Desktop/pats/edfs/PUPbeta/External/');
original = '/Users/sq566/Desktop/pats/inspect/pats-814226-baseline.edf';

[header, signalHeader, signalCell] = blockEdfLoad(original);

plot(signalCell{1,4})
