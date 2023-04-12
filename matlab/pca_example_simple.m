% a clear MATLAB workspace is a clear mental workspace
close all; clear, clc

%PCA on simulated data
loop_thr = linspace(0.1,1,10);

for idx=1:length(loop_thr)
    % data
    x = [ 1*randn(1000,1) loop_thr(idx)*randn(1000,1) ];
    
    % rotation matrix
    th = pi/4;
    R1 = [ cos(th) -sin(th); sin(th) cos(th) ];
    
    % rotate data to induce correlation
    y = x*R1;
    
    % plot the data in a scatter plot
    figure(1), clf
    set(gcf,'Position',[300 300 1800 400]) % set(0,'DefaultFigureWindowStyle','docked')
    % subplot(243)
    plot(y(:,1),y(:,2),'m.','markersize',17)
    
    % make the plot look nicer
    set(gca,'xlim',[-1 1]*max(y(:)),'ylim',[-1 1]*max(y(:)))
    xlabel('C3'), ylabel('C4')
    axis square
    
    %-------------- PCA via eigendecomposition -------------------------
    
    % mean-center ( Data must have zero mean before computing covariance )
    y = y - mean(y,1);
    
    % covariance matrix
    covmat = y'*y / (length(y)-1);
    covmat = round(covmat, 3);
    
    % eigendecomposition
    [evecs,evals] = eig(covmat);
    [evals,sidx] = sort(diag(evals),'descend');
    evecs = evecs(:,sidx);
    
    % plot the eigenvectors on the data
    subplot(131)
    plot(y(:,1),y(:,2),'m.','markersize',17)
    
    % make the plot look nicer
    set(gca,'xlim',[-1 1]*max(y(:)),'ylim',[-1 1]*max(y(:)))
    xlabel('C3'), ylabel('C4')
    axis square
    hold on
    plot([0 evecs(1,1)]*evals(1),[0 evecs(2,1)]*evals(1),'k','linew',3)
    plot([0 evecs(1,2)]*evals(end),[0 evecs(2,2)]*evals(end),'k','linew',3)
    
    % Display eigen vectors
    % ax1 = subplot(245);
    %xCord = repmat(1:2,2,1);
    %yCord = xCord';
    %t = num2cell(evecs);
    %imagesc(evecs), axis square
    %text(xCord(:), yCord(:), t, 'HorizontalAlignment', 'Center')
    %set(gca,'YTickLabel',[], 'XTickLabel',[])
    %title('Spatial filter (Eigen vectors)', 'FontSize',12)
    
    % Component time series
    xCord = repmat(1:2,2,1);  % generate x-coordinates
    yCord = xCord';           % generate y-coordinates
    % subplot(131);

    % Projection
   
    compTS = y * evecs;       % Dot product (If the degree angle between two vector is less than 90, the resulting
                              % dot product will be positive.The dot product is a single 
                              % number that provides information about the relationship between two vectors)   
    PC1 = compTS(:,1);
    PC2 = compTS(:,2);
    
    newY = y;
    newY_covmat = newY'*newY / (length(newY)-1);
    newY_covmat = round(newY_covmat, 3);
    t = num2cell(newY_covmat);     % extact values into cells
    %imagesc(newY_covmat), axis square
    %colormap(summer)
    %text(xCord(:), yCord(:), t, 'HorizontalAlignment', 'Center')
    %set(gca,'clim',[-1 1],'xtick',1:4,'ytick',1:4, 'XAxisLocation', 'top')
    %xticklabels({'C3','C4'})
    %yticklabels({'C3','C4'})
    %xlabel('Channels'), ylabel('Channels')
    %title('Covariance Matrix', 'FontSize',12)
    %clear xCord yCord;
    
    
    % Plot signal
    % subplot(241)
    %plot(y(:,1), 'Color', "#A2142F")
    %title('C3')
    %ylim([-4 4])
    % subplot(242)
    %plot(y(:,2), 'Color', "#77AC30")
    %title('C4')
    %ylim([-4 4])
    subplot(133)
    p1 = plot(PC1, 'Color', "#EDB120");
    hold on
    p2 = plot(PC2, 'Color', "#0072BD");
    p1.Color(4) = 0.75;
    p2.Color(4) = 0.25;
    legend('PC1', 'PC2')
    title('Component time series')
    ylim([-4 4])
    
    % Transform eigen values to an interpretable scale
    evals = 100*evals/sum(evals);
    
    % Plot Eigen spectrum
    
    subplot(132)
    p3 = bar(evals,'FaceColor',[1 0 0],'EdgeColor',[0 0 0],'LineWidth',1.5);
    p3.FaceAlpha = 0.25;
    xticks([1 2])
    xticklabels({'PC1', 'PC2'})
    ylabel('Percent Variance')
    title('Eigen Spectrum')
    text([1,2],evals(:),compose('%f',evals(:)),'HorizontalAlignment','center','VerticalAlignment','bottom')
    ylim([min(ylim),min(ylim)+range(ylim)*1.05])
    fname = ['/Users/sq566/Desktop/' num2str(idx) '.png'];
    print(figure(1),fname,'-r300','-dpng');
end

% --------------- Important notes --------------------
% 1) var(y*evecs) is the same as evecs'*covmat*evecs;
% 2) covmat*evecs is the same as evecs*evals
%    Intuition behind multiplying covmat with evecs
%    The eigen vectors in evecs points to the biggest direction of variance
%    in covmat (covariance matrix)

% ------------------- Plot 3d ------------------------
% figure(2)
% a1 = y(:,1);
% b1 = y(:,2);
% c1 = PC1;
% cmp = linspace(min(a1),max(b1),numel(c1));
% scatter3(a1,b1,c1, [], cmp,'filled'); 
% colorbar;
% view(-80,15)
% xlabel('C3'); ylabel('C4'); zlabel('PC1')
% grid on;grid minor 

