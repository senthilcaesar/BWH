function positioncodes = PositionCodeLookup(poscodes,protocol)

        for i = 1:size(poscodes,1)
            if strcmp(protocol,poscodes{i,1})
                positioncodes = cell2mat(poscodes(i,2:end));%[poscodes{i,2} poscodes{i,3} poscodes{i,4} poscodes{i,5} poscodes{i,6}];
            end
        end
        