function peerfolders=ea_getmagetpeers(options)

switch options.prefs.machine.normsettings.maget_peerset
    case 'Peers from selected cohort'
        if length(options.uipatdirs)<3
            ea_error('Please select more than 1 subject/patient when selecting "Peers from selected cohort" in the MAGeT Normalization settings.');
        end
        peerfolders=options.uipatdirs;
    case 'IXI-Dataset'
        ptage=ea_getpatientage([options.root,options.patientname,filesep]);
        peerfolders=ea_getIXI_IDs(21,ptage);
    case 'User-Defined'
        peerfolders=options.normalize.settings.peersetcell;
end

if length(peerfolders)<3
    ea_error('Please specify a valid number of peers for MAGeT-Brain like approaches');
end
% make sure no peer is also the subject:
ps=ismember(peerfolders,[options.root,options.patientname]);
peerfolders(ps)=[];