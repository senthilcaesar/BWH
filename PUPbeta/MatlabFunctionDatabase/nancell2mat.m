function y = nancell2mat(x)
    
    y = nan(size(x));
    for j=1:size(x,2)
        for i=1:size(x,1)
            try
                temp = x{i,j};
                if isnumeric(temp)
                    y(i,j)=temp; 
                end
            catch
            end
        end
    end
    
    %code by SS, couldn't find simpler way of handling cell2mat errors with
    %text