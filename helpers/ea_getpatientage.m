function ptage=ea_getpatientage(directory)
% this function assumes a text or ml file called ea_age.txt or
% ea_age.mat stating the age of the patient within directory.


if exist([directory,'ea_age.txt'],'file')
    ptage=load([directory,'ea_age.txt']);
    return
end
if exist([directory,'ea_age.mat'],'file')
    ptage=load([directory,'ea_age.mat']);
    fn=fieldnames(ptage);
    ptage=ptage.(fn{1});
return
end

% if non exists return an average age defined in ea_prefs
prefs=ea_prefs('');
ptage=prefs.ixi.meanage;