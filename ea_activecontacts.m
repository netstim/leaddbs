function S = ea_activecontacts(S)

for pt=1:length(S)
    for side=1:2
        switch side
            case 1
                sidec='R';
            case 2
                sidec='L';
        end

        S(pt).activecontacts{side} = zeros(1,S(pt).numContacts);

        for source=1:4
            for c=1:S(pt).numContacts
                try    % struct may not even be specified yet. In this case just keep marked as inactive.
                    if S(pt).([sidec,'s',num2str(source)]).(['k',num2str(c)]).perc>0 && S(pt).([sidec,'s',num2str(source)]).amp>0
                        S(pt).activecontacts{side}(c)=1;
                    end
                end
            end
        end
    end
end
