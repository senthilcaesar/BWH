function PosSignalTranslated = PositionTranslator(positioncodes,positioncodesout,PosSignal)

    PosSignal=PosSignal(:);
    positioncodes=positioncodes(:)';
    error = PosSignal-positioncodes;
    [~,i]=min(abs(error'));
    PosSignalTranslated = positioncodesout(i);
    PosSignalTranslated=PosSignalTranslated(:);
    PosSignalTranslated(isnan(PosSignal))=5;
    