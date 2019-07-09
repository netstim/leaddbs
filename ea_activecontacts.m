function S = ea_activecontacts(S)
% 
%
% USAGE:
%
%    S = ea_activecontacts(S)
%
% INPUT:
%    S:     
%
% OUTPUTS:
%    S:     
%
% .. AUTHOR:
%       - Andreas Horn, Original file
%       - Ning Fey, Original file
%       - Daniel Duarte, Documentation

for pt=1:length(S)
    for side=1:2
        switch side
            case 1
                sidec='R';
                ks=0:7;
            case 2
                sidec='L';
                ks=8:15;
        end
        cs(ks+1)=1:8;

        S(pt).activecontacts{side}=zeros(1,8);

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
