function activeidx = ea_activeidx(S,side,conts,elspec)

% Load contact locations and assign active contact
if side == 1
    sidec = 'R';
else
    sidec = 'L';
end

for con = 1:S.numel
    for source=1:4
        if S.([sidec,'s',num2str(source)]).amp % then this active contact could be from this source since source is active
            if S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).perc % current captured contact is from this source
                activeidx(source).con(con).ix=conts.contactidx{con};
                activeidx(source).con(con).pol=S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).pol;
                activeidx(source).con(con).perc=S.([sidec,'s',num2str(source)]).(['k',num2str(con)]).perc;
            end
        end
    end
end
