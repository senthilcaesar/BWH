%%
%     COURSE: PCA and multivariate neural signal processing
%    SECTION: Dimension reduction with PCA
%      VIDEO: PCA intuition with scatter plots and covariance surfaces
% Instructor: sincxpress.com
%
%%

% a clear MATLAB workspace is a clear mental workspace
close all; clear, clc

%% PCA on simulated data

% data
x = [ 1*randn(1000,1) 1*randn(1000,1) ];

% rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];

% rotate data to induce correlation
y = x*R1;


% plot the data in a scatter plot
figure(1), clf
set(gcf,'Position',[300 300 1800 800])
subplot(231)
plot(y(:,1),y(:,2),'m.','markersize',17)

% make the plot look nicer
set(gca,'xlim',[-1 1]*max(y(:)),'ylim',[-1 1]*max(y(:)))
xlabel('C3'), ylabel('C4')
axis square

%% PCA via eigendecomposition

% mean-center
y = y - mean(y,1);

% covariance matrix
covmat = y'*y / (length(y)-1);

% eigendecomposition
[evecs,evals] = eig(covmat);
[~,maxcomp] = sort(diag(evals));


% plot the eigenvectors on the data
hold on
plot([0 evecs(1,1)]*evals(1),[0 evecs(2,1)]*evals(1),'k','linew',3)
plot([0 evecs(1,2)]*evals(end),[0 evecs(2,2)]*evals(end),'k','linew',3)


% show the covariance matrix
xCord = repmat(1:2,2,1);  % generate x-coordinates
yCord = xCord';           % generate y-coordinates
t = num2cell(covmat);     % extact values into cells
subplot(232)
imagesc(covmat), axis square
text(xCord(:), yCord(:), t, 'HorizontalAlignment', 'Center')
set(gca,'clim',[-1 1],'xtick',1:2,'ytick',1:2, 'XAxisLocation', 'top')
xticklabels({'C3','C4'})
yticklabels({'C3','C4'})
xlabel('Channels'), ylabel('Channels')
title('Covariance Matrix', 'FontSize',12)
clear xCord yCord t;

% Component time series
xCord = repmat(1:4,4,1);  % generate x-coordinates
yCord = xCord';           % generate y-coordinates
subplot(233)
comp1 = y * evecs(:,1);
comp2 = y * evecs(:,maxcomp(end));
newY = y;
newY(:,3) = comp1;
newY(:,4) = comp2;
newY_covmat = newY'*newY / (length(newY)-1);
t = num2cell(newY_covmat);     % extact values into cells
imagesc(newY_covmat), axis square
text(xCord(:), yCord(:), t, 'HorizontalAlignment', 'Center')
set(gca,'clim',[-1 1],'xtick',1:4,'ytick',1:4, 'XAxisLocation', 'top')
xticklabels({'C3','C4', 'Comp1', 'Comp2'})
yticklabels({'C3','C4', 'Comp1', 'Comp2'})
xlabel('Channels'), ylabel('Channels')
title('Covariance Matrix', 'FontSize',12)
colorbar()


% Plot signal
subplot(234)
plot(y(:,1), 'Color', "#A2142F")
title('C3')
ylim([-4 4])
subplot(235)
plot(y(:,2), 'Color', "#77AC30")
title('C4')
ylim([-4 4])
subplot(236)
p1 = plot(comp1, 'Color', "#EDB120");
hold on
p2 = plot(comp2, 'Color', "#0072BD");
p1.Color(4) = 0.75;
p2.Color(4) = 0.25;
legend('Comp1', 'Comp2')
title('Component time series')
ylim([-4 4])
fname = '/Users/sq566/Desktop/10.png';
print(figure(1),fname,'-r600','-dpng');


% Transform eigen values to an interpretable scale
evals = 100*evals/sum(evals);
