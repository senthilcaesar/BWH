function ArSig = ArSignalEffects(ArSig,Opt,Parameter)

% ArSig = 1*(Ensembles.ArPr*1 > 0.5);
% AR1 = ArSig;

if isempty(Parameter)
    Parameter=10;
end

if Opt>0
    AR3=ArSig(:);
    try
        if sum(AR3)>0
            
            Iar = find([NaN;diff(AR3)]==1); %index of arousal onset
            Iar2 = find(diff(AR3)==-1); %index of pre arousal onset
            if Iar2(1)<Iar(1) %tidy up in case they don't match
                Iar = [1;Iar];
            end
            if Iar2(end)<Iar(end)
                Iar2 = [Iar2;length(AR3)];
            end
            lengtharousalB=Iar2-Iar+1;
            for i=1:length(Iar)
                if Opt==2 %exponential decay
                    
                    tdata = (0:lengtharousalB(i)-1);
                    VRAtau=Parameter;
                    RefDur=15;
                    data=(exp(-1/VRAtau*tdata)-exp(-1/VRAtau*lengtharousalB(i)))/(1-exp(-1/VRAtau*RefDur));
                elseif Opt==1 %linear
                    data=(lengtharousalB(i):-1:1)/lengtharousalB(i); %linear reduction with no effect at offset
                end
                AR3(Iar(i):(Iar(i)+lengtharousalB(i)-1))=data; %save signal as appropriate in AR3
            end
        end
        ArSig=AR3;
    catch me
    end
end

