function [Position,PositionRaw] = getPos(SigT,ChannelsList,settings)
PosChan = find(strcmp(ChannelsList,'Position')==1);
if ~isempty(PosChan)
    PositionRaw=SigT.Position;
    Position = PositionTranslator(settings.positioncodes,settings.positioncodesout,PositionRaw);
else  % DLM added this case to handle NaN ColumnHeads(10),
    % It artificially sets the entire night to supine
    if length(settings.supinepositioncode) > 1
        supinecode = settings.supinepositioncode(2);
    else
        supinecode = settings.supinepositioncode;
    end
    Position=ones(height(SigT),1)*supinecode;
    PositionRaw = Position; %to avoid code to break if no position channel available
end