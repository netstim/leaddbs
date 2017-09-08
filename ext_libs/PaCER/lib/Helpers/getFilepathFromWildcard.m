%% getFilepathFromWildcard
%
% Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
function filepath=getFilepathFromWildcard(folder, wildcard)
%     theFile=[folder filesep wildcard];
%     filepath=dir(theFile);
    wildcardRegexp = strrep(wildcard, '*', '.*');
    allFilepaths=dir(folder);
    matchingIndices = [];
    for i=1:numel(allFilepaths)
        fpName = allFilepaths(i).name;
        currIdx = regexpi(fpName, wildcardRegexp);
        if ( ~isempty(currIdx) )
            matchingIndices(end+1) = i;
        end
    end
    filepath = [];
    if ( numel(matchingIndices) >= 1 )
        for i=1:numel(matchingIndices)
            currFp = allFilepaths(matchingIndices(i));
            if ( ~strcmp(currFp.name(1) ,'.'))
                filepath = [folder filesep currFp.name];  
                if ( size(filepath,1) > 1 )
                    warning(['Multiple matches in folder ' folder ...
                        ' for wildcard ' wildcard ...
                        ' ! The first match that does not begin with . is returned ( '...
                        filepath ' ).']);
                end
                return;
            end
        end
    end
%     if ( ~isempty(filepath) )
%         
%         for i=1:numel(filepath)
%             if ( ~strcmp(filepath(i).name(1) ,'.') )
%                 filepath = [folder filesep filepath(i).name];  
%                 if ( size(filepath,1) > 1 )
%                     warning(['Multiple matches in folder ' folder ...
%                         ' for wildcard ' wildcard ...
%                         ' ! The first match that does not begin with . is returned ( '...
%                         filepath ' ).']);
%                 end
%                 return;
%             end
%         end
%     else
%         filepath = [];
%     end
end