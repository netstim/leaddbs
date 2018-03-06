function [changedstates,ret]=ea_checkfschanges(resultfig,fibersfile,seedfile,targetsfile,thresh,mode)
% small helper function that determines changes in fibertracking results.
ret=1;
ofibersfile=getappdata(resultfig,'fibersfile'); % mode independent
oseedfile=getappdata(resultfig,[mode,'seedfile']);
otargetsfile=getappdata(resultfig,[mode,'targetsfile']);
othresh=getappdata(resultfig,[mode,'thresh']);

changedstates=[~isequal(fibersfile,ofibersfile)
    ~isequal(seedfile,oseedfile)
    ~isequal(targetsfile,otargetsfile)
    ~isequal(thresh,othresh)];
if any(changedstates)
    ret=0;
end

setappdata(resultfig,'fibersfile',fibersfile); % mode independent
setappdata(resultfig,[mode,'seedfile'],seedfile);
setappdata(resultfig,[mode,'targetsfile'],targetsfile);
setappdata(resultfig,[mode,'thresh'],thresh);
