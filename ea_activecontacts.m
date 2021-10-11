function S = ea_activecontacts(S)

for pt=1:length(S)
    for side=1:2
        switch side
            case 1
                sidec='R';
                % Aleva
                if length(fieldnames(S.Rs1))==15
                    ks=0:11;
                else
                    ks=0:7;
                end
            case 2
                sidec='L';
                % Aleva
                if length(fieldnames(S.Ls1))==15
                    ks=12:23;
                else
                    ks=8:15;
                end
        end
        % Aleva
        if length(fieldnames(S.Rs1))==15
           cs(ks+1)=1:12;
           S(pt).activecontacts{side}=zeros(1,12);
        else
           cs(ks+1)=1:8;
           S(pt).activecontacts{side}=zeros(1,8);
        end

        for source=1:4
            for cnt=ks
                try    % struct may not even be specified yet. In this case just keep marked as inactive.
                    if S(pt).([sidec,'s',num2str(source)]).(['k',num2str(cnt)]).perc>0 && S(pt).([sidec,'s',num2str(source)]).amp>0
                        S(pt).activecontacts{side}(cs(cnt+1))=1;
                    end
                end
            end
        end
    end
end
