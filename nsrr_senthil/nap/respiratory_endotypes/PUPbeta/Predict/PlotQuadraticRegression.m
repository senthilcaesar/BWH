% 
% PlotSliceOption=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%% Laura: I like 4, 5, 6, 7(?), 11, and 3.3 (3D) as informative 
% switch PlotSliceOption
%     case 1, IthLinear = [1 2]; %xlims = [0 103]; ylims = [0.21 0.78];
%     case 2, Ilist = [1 3 4]; %xlims = [0 103]; ylims = [0.21 0.78];
%     case 3, Ilist = [1 2 4]; %xlims = [0 103]; ylims = [0.21 0.78]; %%%%%
%     case 3.2, Ilist = [1 2 5]; %xlims = [0 103]; ylims = [0.21 0.78]; %%%%%
%     case 3.3, Ilist = [1 4 5]; %xlims = [0 103]; ylims = [0.21 0.78]; %%%%%
%     case 4, Ilist = [1 2]; %xlims = [0.21 0.78]; ylims = [-100 38];
%     case 5, Ilist = [1 4]; %xlims = [0 103]; ylims = [90 200];
%     case 6, Ilist = [2 4]; %xlims = [0.21 0.78]; ylims = [90 200];
%     case 7, Ilist = [1 3]; %xlims = [0 103]; ylims = [-100 38];
%     case 8, Ilist = [2 5]; %xlims = [0.21 0.78]; ylims = [0 100];
%     case 9, Ilist = [1 5]; %xlims = [0.21 0.78]; ylims = [0 100];
%     case 10, Ilist = [3 5]; %xlims = [0.21 0.78]; ylims = [0 100];
%     case 11, Ilist = [4 5]; %xlims = [0.21 0.78]; ylims = [0 100];
%     case 21, Ilist = [1 2 6]; %xlims = [0 103]; ylims = [0.21 0.78];
%     case 22, Ilist = [2 3 6]; %xlims = [0 103]; ylims = [0.21 0.78];
%     case 23, Ilist = [1 3 6]; %xlims = [0 103]; ylims = [0.21 0.78]; %%%%% 
%     case 24,  Ilist = [2 4 5];
% end
if ~exist('markersiz')
markersiz=50;
    end


M=length(xvalueslist);
x=Ain(:,1:M);
xunscaled=Amatrix(:,1:M);
xunscaled(Exclude==1,:)=NaN;
B=Btemp(2:end);


%requires standardize data; n.b. do not re-standardize after generating quadratic terms. 

Is = zeros(M,1);
Is(Ilist)=1;

K=Btemp(1); %model constant

%linear terms
L = zeros(M,1);
for i=1:length(L)
    I=find(Ilisttest==i);
    if ~isempty(I)
        L(i)=B(I);
    end
end
L=L(Ilist);

Qlong = zeros((M^2-M)/2+M,1);
for i=1:length(Qlong)
    I=find(Ilisttest==i+M);
    if ~isempty(I)
        Qlong(i)=B(I);
    end
end

Q=0*Q0;
count=1;
for j=1:M
    for i=1:M
        if Q0(i,j)~=0
            Q(i,j)=Qlong(count);
            count=count+1;
        end
    end
end
if 0 %seems unnecessary, already are zeros
    Q(linearonlyterms,:)=0;
    Q(:,linearonlyterms)=0;
end
Q=Q(Ilist,Ilist);


myf2 = @(A) K + [A]*L + sum(([A]*Q') .* [A], 2);
  
PredtTemp = myf2(x(:,Ilist));

[PredtTemp PrPredict];
%figure(5); plot(PredtTemp-thresopt,PrPredict-thresopt,'.')

%thresopt=25

PredTsimp = 1*(PredtTemp>thresopt);
    PredTsimp(Exclude==1)=NaN;

predmodel.model=myf2;
predmodel.K=K-thresopt;
predmodel.L=L*2;
predmodel.Q=Q;

mismatch=(PredTsimp~=(PrPredict>thresopt));
        %Imismatch=[];
%mismatch=0*mismatch
plotrange=[1:size(Amatrix,1)]; %plot Cistulli + Marques 
I=(Exclude==0); sum(I);
plotrange(Exclude==1)=[];
plotrangeY=Exclude==0&mismatch==0;

%
C(plotrangeY==0,:)=0.7;

%

if length(Ilist)==3
    
    
%figure(4); clf(4);


xcolor = 1;
            %cutoff_=50;
            %criteriaRplot = Yvariable>cutoff_;
            %criteriaN = Yvariable<=cutoff_;     
        
        %NaN outside plotrange
        
        
            % 3d figure in spheres
            compressxyz=[1 1 1]; %[1 0.5 0.5]
            
            set(gcf,'color',[1 1 1]);
            %htemp=scatter3(Amatrix(:,Ilist(1)),Amatrix(:,Ilist(2)),Amatrix(:,Ilist(3)));
            htemp=scatter3(Amatrix(:,Ilist(1)),Amatrix(:,Ilist(2)),Amatrix(:,Ilist(3)),'marker','none');
            xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
            yrange = [min(Amatrix(plotrange,Ilist(2))) max(Amatrix(plotrange,Ilist(2)))];
            zrange = [min(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3)))];
            cubevol = diff(xrange)*diff(yrange)*diff(zrange);
            %ballsize=0.02*cubevol;
            ballsize=1*max([diff(xrange) diff(yrange) diff(zrange)]);
            %ballsize=0.02*cubevol^(1/3);
            camerapos=get(gca,'CameraPosition');
            
%             criteriaR = Yvariable>66.67;
%             criteriaMid = Yvariable>50&~criteriaR;
%             criteriaMid2 = Yvariable>33.33&~criteriaMid&~criteriaR;
%             criteriaN = (~criteriaR)&(~criteriaMid)&(~criteriaMid2)&~isnan(Yvariable);

%             for i=1:length(Yvariable)
%                 distancefromeye(i)=sqrt((Amatrix(i,Ilist(1))-camerapos(1)).^2+(Amatrix(i,Ilist(2))-camerapos(2)).^2+(Amatrix(i,Ilist(3))-camerapos(3)).^2);
%             end
%             distancefromeye=distancefromeye/nanmean(distancefromeye);
%             distancefromeye=distancefromeye.^1.5;
%             distancefromeye=1./distancefromeye;
            switch offscreenpredoption
                case 1
                    %do nothing
                case 2
                    plotrange=[1:size(Amatrix,1)]; %plot Cistulli + Marques
                    I=(Exclude==0); sum(I)
                    plotrange(Exclude==1)=[];
                    C(mismatch,:)=0.7;
            end

            distancefromeye= 0*Amatrix(:,Ilist(3))+1;
            hold('on');
            h=scatter3sph(Amatrix(plotrange,Ilist(1)),Amatrix(plotrange,Ilist(2)),Amatrix(plotrange,Ilist(3)),'size',ballsize*distancefromeye(plotrange),'color',C(plotrange,:),'transp',1,'compress',compressxyz);
            hold('on');
            xrange = [min(Amatrix(plotrange,Ilist(1))) max(Amatrix(plotrange,Ilist(1)))];
            yrange = [min(Amatrix(plotrange,Ilist(2))) max(Amatrix(plotrange,Ilist(2)))];
            zrange = [min(Amatrix(plotrange,Ilist(3))) max(Amatrix(plotrange,Ilist(3)))];
            
            Fextraxyz=0.02; %extra "white space" variable.
            xlim([min(xrange)-Fextraxyz*diff(xrange) max(xrange)+Fextraxyz*diff(xrange)]);
            ylim([min(yrange)-Fextraxyz*diff(yrange) max(yrange)+Fextraxyz*diff(yrange)]);
            zlim([min(zrange)-Fextraxyz*diff(zrange) max(zrange)+Fextraxyz*diff(zrange)]);
            %surface "wall"
             %lighting for depth
            set(h,'facelighting','phong','ambientstrength',0.5);
            h2=light('position',[nanmean(Amatrix(plotrange,Ilist(1))) nanmean(Amatrix(plotrange,Ilist(2))) nanmean(Amatrix(plotrange,Ilist(3)))],'style','local');
            set(gca,'CameraPosition',camerapos,'projection','perspective','gridlinestyle','-','TickLength',[0 0],'fontsize',14);
            set(gca,'fontname','arial narrow');
            
            xlabel(xvalueslist{Ilist(1)},'interpreter','none');
            ylabel(xvalueslist{Ilist(2)},'interpreter','none');
            zlabel(xvalueslist{Ilist(3)},'interpreter','none');
            
Plot3Dclassifier(predmodel,xunscaled(:,Ilist),shift(Ilist),criteriaR,"Quadratic")

elseif length(Ilist)==2
   %figure(5); clf(5);
   
    trans=0.85;
       %     scatter(xunscaled(plotrangeY==0,Ilist(1)),xunscaled(plotrangeY==0,Ilist(2)),35,0.6*[1 1 1],'markerfacealpha',trans);
            hold('on');
%             scatter(xunscaled(plotrangeY==1&criteriaRplot==1,Ilist(1)),xunscaled(plotrangeY==1&criteriaRplot==1,Ilist(2)),markersiz,[0.9 0.4 0],'filled','markerfacealpha',trans);
%             
%             scatter(xunscaled(plotrangeY==1&criteriaRplot2==1,Ilist(1)),xunscaled(plotrangeY==1&criteriaRplot2==1,Ilist(2)),markersiz,[0 0.6 0],'filled','markerfacealpha',trans);
%             
           % set(gca,'box','off','fontname','arial narrow','fontsize',14);
            
%             scatter(xunscaled(plotrangeY==1&criteriaN2==1,Ilist(1)),xunscaled(plotrangeY==1&criteriaN2==1,Ilist(2)),markersiz,[0.7 0 0],'filled','markerfacealpha',trans); %more clear NR
%             scatter(xunscaled(plotrangeY==1&criteriaN==1,Ilist(1)),xunscaled(plotrangeY==1&criteriaN==1,Ilist(2)),markersiz,[0.7 0 0],'filled','markerfacealpha',trans);
%             
         %   scatter(xunscaled(plotrangeY==1,Ilist(1)),xunscaled(plotrangeY==1,Ilist(2)),markersiz,C(plotrangeY==1,:),'filled','markerfacealpha',trans);
            
            
            Fextra=0.3;
            xrange = [min(xunscaled(plotrange,Ilist(1))) max(xunscaled(plotrange,Ilist(1)))];
                xrange = [xrange(1)-diff(xrange)*Fextra xrange(2)+diff(xrange)*Fextra];
            yrange = [min(xunscaled(plotrange,Ilist(2))) max(xunscaled(plotrange,Ilist(2)))];
                yrange = [yrange(1)-diff(yrange)*Fextra yrange(2)+diff(yrange)*Fextra];
                lims = [xrange yrange];
                [X,Y] = meshgrid(linspace(xrange(1),xrange(2)),linspace(yrange(1),yrange(2)));
                Xorig = X; Yorig = Y;
                
                scaleFactors = 1./nanstd(xunscaled);
                shift = -nanmean(xunscaled,1);
                if 1
                    X = scaleFactors(Ilist(1)) .* (X + shift(Ilist(1)) );
                    Y = scaleFactors(Ilist(2)) .* (Y + shift(Ilist(2)) );
                end
                %svm:
                %Z = (feval(SVMModel.KernelFunction,SVMModel.SupportVectors,[X(:),Y(:)],SVMModel.KernelFunctionArgs{:})'*SVMModel.Alpha(:)) + SVMModel.Bias; 
                Z=myf2([X(:),Y(:)]);
                hold('on');
                %scatter(Xorig(:),Yorig(:),20,[Z(:)>thresopt 1+0*Z(:) 0*Z(:)])
              %  contour(Xorig,Yorig,reshape(Z-thresopt,size(X)),[0 0],'k');
              if 1 % if red background color for pred non-responder - added SO 03/04/2022
                  set(gca,'color',[ 0.95 0.78 0.78])   
              end
              
                contourf(Xorig,Yorig,reshape(Z-thresopt,size(X)),[0 0],'k')
              
               
                %%% set colormap colors so plots green region
                caxis([-1.5 1.1])
                map=      [ 0.4706    0.6706    0.1882
                    0.5765    0.7365    0.3506
                    0.6824    0.8024    0.5129
                    0.7882    0.8682    0.6753
                    0.8941    0.9341    0.8376
                    1.0000    1.0000    1.0000];
                colormap(map); 
               
                
                 scatter(xunscaled(plotrangeY==0,Ilist(1)),xunscaled(plotrangeY==0,Ilist(2)),35,0.6*[1 1 1],'markerfacealpha',trans);
           
            scatter(xunscaled(plotrangeY==1,Ilist(1)),xunscaled(plotrangeY==1,Ilist(2)),markersiz,C(plotrangeY==1,:),'filled','markerfacealpha',trans);
            
           set(gca,'box','off','fontname','arial narrow','fontsize',14);
                
                
                xlabel(xvalueslist(Ilist(1)),'interpreter','none','fontsize',9);
                ylabel(xvalueslist(Ilist(2)),'interpreter','none','fontsize',9);
                
else
    disp('"Plot3D" needs 2 or 3 parameters')
end

%
axisorder={'x','y','z'};
for i=1:length(Ilist)

        if string(xvalueslist(Ilist(i)))==string('VpassiveT')
            VPvals = [0 50 80 95 100];
        VPvals2 = (1-([1-VPvals/100].^0.5))*100;
        set(gca,[axisorder{i} 'tick'],VPvals2,[axisorder{i} 'ticklabels'],VPvals);
        temp='Vpassive, %';
        eval([axisorder{i} 'label(temp)']); 
        end
        
        if string(xvalueslist(Ilist(i)))==string('ArThresT')
        ATvals = [100 110 150 200 250];
        ATvals2 = (1+([ATvals/100-1].^0.5))*100;
        set(gca,[axisorder{i} 'tick'],ATvals2,[axisorder{i} 'ticklabels'],ATvals);
        temp='Arousal Threshold, %';
        eval([axisorder{i} 'label(temp)']);
        eval([axisorder{i} 'lim([80 240])']);
        end
        
        if string(xvalueslist(Ilist(i)))==string('AHIbaselineT')
        AHIvals = [0 5 15 30 60 90];
        AHIvals2 = AHIvals.^0.5;
        set(gca,[axisorder{i} 'tick'],AHIvals2,[axisorder{i} 'ticklabels'],AHIvals);
        temp='AHI, events/hr';
        eval([axisorder{i} 'label(temp)']);
        eval([axisorder{i} 'lim([80 240])']);
        end
end
 set(gca,'tickdir','out')
 set(gcf,'position',[450   350   560   470]);

xlim([xrange(1)-diff(xrange)*0.005 xrange(2)+diff(xrange)*0.005]);
ylim([yrange(1)-diff(yrange)*0.005 yrange(2)+diff(yrange)*0.005]);
%ylim(ylims);

