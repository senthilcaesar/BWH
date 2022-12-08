%add to regression
quadraticterm=nan(1,length(labels));
for iii=1:length(labels)
    quadraticterm(iii)=contains(labels(iii),':')==1 | contains(labels(iii),'^')==1;
end

LinearTerms = ~quadraticterm;


%setup colors  
if ~exist('criteriaPlot')
criteriaPlot = [...
    T.Responder==1 ...
    isnan(T.Responder) ...
    T.Responder==0 ...
    T.Responder==0 ...
    ];
end

if ~exist('colorscheme')
    colorscheme=5;
end

switch colorscheme
    case 1 %dark green, light green, orange, red
        Coptions = [0 0.5 0;0.35 0.7 0;0.9 0.5 0;0.75 0.1 0.1];
    case 2 %dark green, yellow, orange, red
        Coptions = [0 0.5 0;0.9 0.5 0;0.9 0.25 0;0.75 0.1 0.1];
    case 3 %dark green, orange, red, dark red
        Coptions = [0 0.5 0;0.9 0.5 0;0.75 0.2 0.1;0.5 0 0];
    case 4 %green (R), green (>50%), red (NR), dark red (NR+worse)
        Coptions = [0 0.5 0;0 0.5 0;0.75 0.2 0.1;0.5 0 0];    
    case 5 %dark green (>67%), light green (>50%, %), red (NR), dark red (NR+worse)
        Coptions = [0 0.5 0;0.35 0.7 0;0.75 0.2 0.1;0.5 0 0];    
    case 6 %red, orange (NaN), green,  green
        %Coptions = [0 0.5 0;0.9 0.5 0;0.75 0.2 0.1;0.75 0.2 0.1];         
        Coptions = [0.75 0.2 0.1;0.9 0.5 0;0 0.5 0;0 0.5 0];         
end


maxlinearforplot = 10;

figure(123); clf(123);
set(gcf,'color',[1 1 1]);
IlisttestLinear = Ilisttest(LinearTerms(Ilisttest));
IlisttestLinear((maxlinearforplot+1):end)=[];
for iii=1:(length(IlisttestLinear)-1)
    for jjj=2:(length(IlisttestLinear))
        if iii>=jjj
            continue
        else
            subplot(length(IlisttestLinear)-1,length(IlisttestLinear)-1,iii+(jjj-1-1)*(length(IlisttestLinear)-1));
            Ilist = [IlisttestLinear(iii) IlisttestLinear(jjj)];
            
            C=ones(size(Amatrix,1),3);
            for i=1:size(Amatrix,1)
                row = find(criteriaPlot(i,:)==1,1,'first');
                C(i,:)=Coptions(row,:);
            end
            offscreenpredoption=1; %grey not disappeared in 3D plots

            run PlotQuadraticRegression
            
            xx = get(gca,'position');
            xx(3:4)= xx(3:4)*1.2;
            set(gca,'position',xx);
            set(gca,'fontsize',8,'fontweight','normal');
            if jjj<length(IlisttestLinear)
                set(gca,'XLabel',[],'XTickLabels',[]);
                
            end
            if iii>1
                set(gca,'YLabel',[],'YTickLabels',[]);
            end
        end
    end
end

