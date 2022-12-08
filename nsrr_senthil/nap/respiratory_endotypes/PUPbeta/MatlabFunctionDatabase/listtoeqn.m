function covarstr = listtoeqn(covarlist)
    if isempty(covarlist)
        covarstr = '';
        return
    end


    covarstr = covarlist{1}; 
    for i=2:length(covarlist), 
        covarstr = [covarstr ' + ' covarlist{i}]; 
    end 
