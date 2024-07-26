% the following code reproduces Figure 3 in 
% https://archive.physionet.org/physiotools/mse/tutorial/tutorial.pdf
% john.malik@duke.edu

clear
close all

trials = 30;

scales = 20;

% store the sample entropy values
t = zeros(scales, 1);

r = 0.15;

m = 2;

n = 3000; % the example uses 30,000

for tau = 1:scales
    
    for j = 1:trials
        
        signal = randn(n, 1);
        
        % [ t(tau, j), A, B ] = multiscaleSampleEntropy_compatible( signal, m, r, tau );
        [ t(tau, j), A, B ] = multiscaleSampleEntropy( signal, m, r, tau );
        
    end
    
    
end

% average over all trials
plot(mean(t, 2))