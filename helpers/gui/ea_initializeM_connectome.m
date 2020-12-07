function M=ea_initializeM_connectome
M=struct;
M.patient.list={};
M.patient.group=[];

M.guid=datestr(datevec(now), 'yyyymmddHHMMSS' ); % each lead groupanalysis needs a unique ID for VTA handling / identification.
M.clinical.vars={};
M.clinical.labels={};
M.ui=struct;
M.ui.listselect=1;
M.ui.clinicallist=1;
M.ui.normregpopup=1;
M.ui.lc.graphmetric=1;
M.ui.lc.normalization=1;
M.ui.lc.smooth=1;
