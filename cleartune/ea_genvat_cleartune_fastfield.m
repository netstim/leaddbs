function Vvate=ea_genvat_cleartune_fastfield(S,side,options,stimname,lgfigure,electrode,elt)

conductivity = options.prefs.machine.vatsettings.fastfield_cb;  % 0.1;

resultfig=getappdata(lgfigure,'resultfig');
elstruct=getappdata(resultfig,'elstruct');
options=getappdata(resultfig,'options');
options.usediffusion=0;
setappdata(resultfig,'elstruct',elstruct);

switch side
    case 1
        sidec='R';
        cnts={'k0','k1','k2','k3','k4','k5','k6','k7'};
    case 2
        sidec='L';
        cnts={'k8','k9','k10','k11','k12','k13','k14','k15'};
end

if ~isfield(S, 'sources')
    S.sources=1:4;
end

Efield=zeros(100,100,100);
for source=S.sources

    stimsource=S.([sidec,'s',num2str(source)]);
    % constvol is 1 for constant voltage and 0 for constant current.
    amp1 = stimsource.amp;
    if amp1>0
        count1=1;
        for cnt=1:length(cnts)
            perc(cnt) = stimsource.(cnts{cnt}).perc;
            if perc(cnt)>0
                Im(count1)=stimsource.(cnts{cnt}).imp;
                count1=count1+1;
            end
        end

        constvol=stimsource.va==1;
        if constvol
            amp_mode = 'V';
            impedance = mean(Im)*1000;
        else
            amp_mode = 'mA';
            impedance = [];
        end

        [Efield_to_add] = ea_get_efield(perc,elt.standard_efield,amp1,conductivity,amp_mode,impedance);
        Efield = Efield+Efield_to_add;
    end
end

electrode_patient = elstruct;

[trans_mat] = get_trans_mat(electrode,electrode_patient,elt.grid_vec,side);

gv=elt.grid_vec;

% create nifti
[~, ~, endian] = computer;
switch endian
    case 'L'
        endian = 0;
    case 'B'
        endian = 1;
end

res=100;
chun1=randperm(res); chun2=randperm(res); chun3=randperm(res);
Vvate.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
Vvate.mat = trans_mat * Vvate.mat;
Vvate.dim=[res,res,res];
Vvate.dt = [16, endian];
Vvate.n=[1 1];
Vvate.descrip='lead dbs - vat';

stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
ea_mkdir(stimDir);
filePrefix = ['sub-', options.subj.subjId,'_sim-'];
fileSuffix = ['_stimset-', S.label];
switch side
    case 1
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-cleartunefastfield_hemi-R', fileSuffix, '.nii'];
    case 2
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-cleartunefastfield_hemi-L', fileSuffix, '.nii'];
end

Vvate.img = Efield; % permute(eeg,[2,1,3]);
%ea_write_nii(Vvate); 