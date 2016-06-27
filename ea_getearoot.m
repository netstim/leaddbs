function root=ea_getearoot

% small function determining the location of the lead-dbs root directory.
if isdeployed
    errordlg(ctfroot);
    root=[ctfroot,'Lead_DBS'];
else
    root=[fileparts(which('lead.m')),filesep];
end