%% SS_area overlay
function addSSArea(InspFlow, InspTime, ShadeColor, LineColor)
[Flow_peaks] = FindPeaksInInspTertiles(InspFlow);
switch size(Flow_peaks,1) %AA_NumOfPeaks(i)
    case 1 % do nothing
        
    case 2 % if AA_NumOfPeaks(i) = 2, simple area
        % find the area between flow signal, and a chord between peaks
        x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(2,1))];
        y = [Flow_peaks(1,2) Flow_peaks(2,2)];
        c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
        Xrange = InspTime(Flow_peaks(1,1):Flow_peaks(2,1));
        ln = length(Xrange);
        chordA = NaN(ln,1);
        for k=1:ln; chordA(k) = (c(2)*(Xrange(k)))+c(1); end
        % find abs area between the chord (between peaks) and flow
        %SS_Area(i)=abs(trapz(Xrange,InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))) ...
        %    - trapz(Xrange,chordA));
        diff_line = chordA' - InspFlow(Flow_peaks(1,1):Flow_peaks(2,1));
        diff_line = max(diff_line,0); % only want the positive values
        %plot(Xrange, chordA, LineColor);
        jbfill(Xrange, chordA, InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))', ShadeColor, LineColor, 0, 0.1);
%         x2 = [Xrange, flipud(Xrange)];
%         inbetween = [InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))' flipud(chordA)];
%         fill(x2, inbetween, 'g');
 
    case 3 % if AA_NumOfPeaks(i) = 3, max ( A-B + B-C ; A-C )
        % AB chord (= between peaks 1 and 2)
        x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(2,1))];
        y = [Flow_peaks(1,2) Flow_peaks(2,2)];
        c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
        XrangeAB = InspTime(Flow_peaks(1,1):Flow_peaks(2,1));
        ln = length(XrangeAB);
        chordAB = NaN(ln,1);
        for k=1:ln; chordAB(k) = (c(2)*(XrangeAB(k)))+c(1); end
        % AB area
        %Area_AB=abs(trapz(XrangeAB,InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))) ...
        %    - trapz(XrangeAB,chordAB));
        diff_lineAB = chordAB' - InspFlow(Flow_peaks(1,1):Flow_peaks(2,1));
        diff_lineAB = max(diff_lineAB,0); % only want the postive values
        Area_AB = trapz(XrangeAB,diff_lineAB);
        % BC chord (= between peaks 2 and 3)
        x = [InspTime(Flow_peaks(2,1)) InspTime(Flow_peaks(3,1))];
        y = [Flow_peaks(2,2) Flow_peaks(3,2)];
        c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
        XrangeBC = InspTime(Flow_peaks(2,1):Flow_peaks(3,1));
        ln = length(XrangeBC);
        chordBC = NaN(ln,1);
        for k=1:ln; chordBC(k) = (c(2)*(XrangeBC(k)))+c(1); end
        % BC area
        %Area_BC=abs(trapz(XrangeBC,InspFlow(Flow_peaks(2,1):Flow_peaks(3,1))) ...
        %    - trapz(XrangeBC,chordBC));
        diff_lineBC = chordBC' - InspFlow(Flow_peaks(2,1):Flow_peaks(3,1));
        diff_lineBC = max(diff_lineBC,0); % only want the postive values
        Area_BC = trapz(XrangeBC,diff_lineBC);
        % AC chord (= between peaks 1 and 3)
        x = [InspTime(Flow_peaks(1,1)) InspTime(Flow_peaks(3,1))];
        y = [Flow_peaks(1,2) Flow_peaks(3,2)];
        c = [[1; 1] x(:)]\y(:);   % Calculate Parameter Vector
        XrangeAC = InspTime(Flow_peaks(1,1):Flow_peaks(3,1));
        ln = length(XrangeAC);
        chordAC = NaN(ln,1);
        for k=1:ln; chordAC(k) = (c(2)*(XrangeAC(k)))+c(1); end
        % AC area
        %Area_AC=abs(trapz(XrangeAC,InspFlow(Flow_peaks(1,1):Flow_peaks(3,1))) ...
        %    - trapz(XrangeAC,chordAC));
        diff_lineAC = chordAC' - InspFlow(Flow_peaks(1,1):Flow_peaks(3,1));
        diff_lineAC = max(diff_lineAC,0); % only want the postive values
        Area_AC = trapz(XrangeAC,diff_lineAC);
         % use the greatest value
        [~, indx] = max([(Area_AB + Area_BC) Area_AC]);
        if indx == 1 % we used chords AB and BC
            %plot(XrangeAB, chordAB, 'r--');
            %plot(XrangeBC, chordBC, 'r--');
            
            jbfill(XrangeAB, chordAB, InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))', ShadeColor, LineColor, 0, 0.1);
            jbfill(XrangeBC, chordBC, InspFlow(Flow_peaks(2,1):Flow_peaks(3,1))', ShadeColor, LineColor, 0, 0.1);
            
%             x2 = [XrangeAB, flipud(XrangeAB)];
%             inbetween = [InspFlow(Flow_peaks(1,1):Flow_peaks(2,1))' flipud(chordAB)];
%             fill(x2, inbetween, 'g');
%             x2 = [XrangeBC, flipud(XrangeBC)];
%             inbetween = [InspFlow(Flow_peaks(2,1):Flow_peaks(3,1))' flipud(chordBC)];
%             fill(x2, inbetween, 'g');
        else % we used chordAC
            %plot(XrangeAB, chordAB, 'g:');
            %plot(XrangeBC, chordBC, 'g:');
            %plot(XrangeAC, chordAC, 'r--');
            jbfill(XrangeAC, chordAC, InspFlow(Flow_peaks(1,1):Flow_peaks(3,1))', ShadeColor, LineColor, 0, 0.1);
%             x2 = [XrangeAC, flipud(XrangeAC)];
%             inbetween = [InspFlow(Flow_peaks(1,1):Flow_peaks(3,1))' flipud(chordAC)];
%             fill(x2, inbetween, 'g');
        end
    otherwise
        % do nothing
end
end