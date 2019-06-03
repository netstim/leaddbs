function [fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_calcdiscfibers(M,discfiberssetting)

I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);
try
    efieldthresh=discfiberssetting.efieldthresh;
catch
    efieldthresh=150;
end
tic

% Get discriminative fiber setting
connthreshold = discfiberssetting.connthreshold/100;
statmetric = discfiberssetting.statmetric;

% protocol selection to be able to check if same analysis has been run
% before.
% opts.percent=connthreshold; % need not protocol anymore, is now dynamic
% option.
opts.patientselection=M.ui.listselect;
opts.regressor=I;
opts.connectome=M.ui.connectomename;
opts.allpatients=M.patient.list;
opts.mirrorsides=M.ui.mirrorsides;
opts.statmetric=statmetric;

if M.ui.mirrorsides
    msuffix='_mirrored';
else
    msuffix='';
end

switch statmetric
    case 1 % ttests
        savesuffix='_ttests';
    case 2 % spearmans R
        savesuffix='_spearmansrho';
end

[reforce,connectomechanged]=ea_discfibers_checkpresence(M,opts); % only static opts need to be equal.
if reforce
    allroilist=cell(length(M.patient.list),2);
    switch statmetric
        case 1 % use paired T-Tests and binary VTA
            suffix='';
        case 2 % use Spearman Rs and E-Fields
            suffix='_efield';
            prefs=ea_prefs;
            if strcmp(prefs.lcm.vatseed,'efield_gauss')
                suffix='_efield_gauss';
            end
    end
    for sub=1:length(M.patient.list) % all patients - for connected fibers selection
        if ~M.ui.mirrorsides
            allroilist{sub,1}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_right.nii'];
            allroilist{sub,2}=[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_left.nii'];
        else
            allroilist(sub,:)=ea_genflippedjointnii([M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_right.nii'],[M.patient.list{sub},filesep,'stimulations',filesep,'gs_',M.guid,filesep,'vat',suffix,'_left.nii'],statmetric==1);
        end
    end
    if ~exist([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'file')
        cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
    else
        if connectomechanged
            cfile=[ea_getconnectomebase('dMRI'),M.ui.connectomename,filesep,'data.mat'];
        else
            cfile=[M.ui.groupdir,'connected_fibers',msuffix,'.mat'];
        end
    end
    switch statmetric
        case 1 % ttests
            [fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_heatfibertracts(cfile,{allroilist},M.ui.listselect,{I},connthreshold);
        case 2 % spearmans R
            [fibsweighted,fibsin,fibsval,iaix]=ea_discfibers_heatfibertracts_corr(cfile,{allroilist},M.ui.listselect,{I},efieldthresh);
    end
    save([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'fibsin','opts','-v7.3');
    save([M.ui.groupdir,'correlative_fibertracts_fibsval',msuffix,savesuffix,'.mat'],'fibsval','iaix','-v7.3');
    save([M.ui.groupdir,'correlative_fibertracts',msuffix,savesuffix,'.mat'],'fibsweighted','opts','-v7.3');
else
    load([M.ui.groupdir,'connected_fibers',msuffix,'.mat'],'fibsin','opts');
    load([M.ui.groupdir,'correlative_fibertracts_fibsval',msuffix,savesuffix,'.mat']);
    load([M.ui.groupdir,'correlative_fibertracts',msuffix,savesuffix,'.mat']);
end
