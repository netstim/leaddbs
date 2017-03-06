function root=ea_getearoot

% small function determining the location of the lead-dbs root directory.
if isdeployed
    root=[ctfroot,filesep];
else
    root=[fileparts(which('lead.m')),filesep];
end