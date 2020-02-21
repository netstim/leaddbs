%% getFilenameFromWildcard
%
% Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
function filepath=getFilenameFromWildcard(folder, wildcard)
    theFile=[folder filesep wildcard];
    filepath=dir(theFile);
    if ( ~isempty(filepath) )
        filepath = [filepath(1).name];                
    else
        filepath = [];
    end
end