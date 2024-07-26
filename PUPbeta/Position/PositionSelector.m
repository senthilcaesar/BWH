function Poscriteria = PositionSelector(PosSignalTranslated,Selection)

%UniversalPosCodes = [Supine	Left	Right	Prone	Unknown	  Upright]

%%
switch Selection
    case 'All'
        disp('"All positions" selected')
        Poscriteria = ones(length(PosSignalTranslated),1);
    case 'NonSupine' % left, right, prone
        codePos = [2 3 4];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;
    case 'Lateral' % left, right
        codePos = [2 3];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;   
    case 'Left'
        codePos = [2];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;
    case 'Right'
        codePos = [3];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;    
    case 'Supine'
        codePos = [1];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;    
    case 'SupineAndNaN'
        codePos = [1 5];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;            
    case 'Prone'
        codePos = [4];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;    
    case 'Flat' % everything excluding upright and unknown
        codePos = [1 2 3 4];
        Poscriteria = sum(PosSignalTranslated == codePos,2)>0;  
    otherwise 
        disp('Position selection error - using all positions')
        Poscriteria = ones(length(PosSignalTranslated),1);
end