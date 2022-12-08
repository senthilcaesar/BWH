% Testing 
%   how do flow and pnasal compare in GS Vs Pred
%   what features contribute to the output?

GS = Gtest_matched;         % 55250, now 63059
Step = 25; % 3 also works pretty well

%PredF = predyL1O_array_flow_matched(:,Step);
%PredP = predyL1O_array_pnasal_matched(:,Step);
PredF = predyL1O_array_flow(:,Step);
PredP = predyL1O_array(:,Step); % in the TnT pnasal data this is pnasal

figure(1); clf(figure(1));
subplot(1,2,1);
scatter(GS, PredF, 2, 'filled','markerfacealpha',0.2);
title('Flow'); xlabel('GS'); ylabel('Pred');
subplot(1,2,2);
scatter(GS, PredP, 2, 'filled','markerfacealpha',0.2);
title('Pnasal'); xlabel('GS'); ylabel('Pred');

for FtrNum=1:108
    
    FtrF = Amatrix2_flow_matched(:,FtrNum); % 63059
    FtrP = Amatrix2_pnasal_matched(:,FtrNum);
    
    figure(3); clf(figure(3));
    subplot(1,2,1);
    scatter(GS, FtrF, 2, 'filled','markerfacealpha',0.2);
    title('Flow'); xlabel('GS'); ylabel('Ftr');
    subplot(1,2,2);
    scatter(GS, FtrP, 2, 'filled','markerfacealpha',0.2);
    title('Pnasal'); xlabel('GS'); ylabel('Ftr');
    
    disp('');
    
end
