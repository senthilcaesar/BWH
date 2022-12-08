figure(1); %phenotype plot
            h=subplot(1,3,1);
            temp=get(h,'children');
            
            color1=[0.9500 0.1000 0.1];
            color2=[0.100 0.4000 0.7];
            color3=[0.500 0.3000 0.4];
            set(temp(6),'FaceColor',color1);
            if settings.runcompare
                set(temp(15),'FaceColor',color2);
            if settings.runcompare2
                set(temp(24),'FaceColor',color3);
            end
            end
            %                 if positioncompare
            %                 set(temp(15),'FaceColor',[0.9500 0.1000 0.1]);
            %                 end
            %
            ranges = [3 4 5 6 9 10 11]; %histograms
            for ii=1:length(ranges)
                h=subplot(2,6,ranges(ii));
                temp=get(h,'children');
                
                %find handles for the histograms only
                clear remove
                for j=1:length(temp)
                    remove(j)=0;
                    if temp(j).Type=="line"
                        remove(j)=1;
                    end
                end
                temp(remove==1)=[];
                
                set(temp(1),'FaceColor',color1);
                if settings.runcompare
                    set(temp(2),'FaceColor',color2);
                if settings.runcompare2
                    set(temp(3),'FaceColor',color3);
                end
                end
                %                 if positioncompare
                %                     set(temp(2),'FaceColor',[0.9500 0.1000 0.1]);
                %                 end
            end
            
            for ii=1:length(ranges)
                h=subplot(2,6,ranges(ii));
                temp=get(h,'children');
                
                %find handles for the histograms only
                clear remove
                for j=1:length(temp)
                    remove(j)=0;
                    if temp(j).Type=="histogram"
                        remove(j)=1;
                    end
                end
                temp(remove==1)=[];
                
                set(temp(1),'Color',color1,'linewidth',3);
                if settings.runcompare
                    set(temp(2),'Color',color2,'linewidth',3);
                if settings.runcompare2
                    set(temp(3),'Color',color3,'linewidth',3);
                end
                end
                %                 if positioncompare
                %                     set(temp(2),'FaceColor',[0.9500 0.1000 0.1]);
                %                 end
            end
            
            newpos = 1.0e+03 * [0.0103    0.2557    1.2613    0.3447];
            set(gcf,'position',newpos)
            
            figure(4); %phenotype plot
            try
                temp=get(gca,'children');
                set(temp(6),'FaceColor',color1);
                set(temp(13),'FaceColor',color1);
                if settings.runcompare
                    set(temp(6+14),'FaceColor',color2);
                    set(temp(13+14),'FaceColor',color2);
                    if settings.runcompare2
                        set(temp(6+28),'FaceColor',color3);
                        set(temp(13+28),'FaceColor',color3);
                    end
                end
            catch
                
            end
            %                 if positioncompare
            %                 set(temp(15),'FaceColor',[0.9500 0.1000 0.1]);
            %                 end