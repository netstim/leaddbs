function varargout=ea_genvat_fastfield(varargin)

useSI = 1;

if nargin==5
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
elseif nargin==6
    acoords=varargin{1};
    S=varargin{2};
    side=varargin{3};
    options=varargin{4};
    stimname=varargin{5};
    lgfigure=varargin{6};
elseif nargin==1
    if ischar(varargin{1}) % return name of method.
        varargout{1}= 'Fastfield (Baniasadi 2020)';
        varargout{2} = true; % Support directed lead
        return
    end
end

conductivity = options.prefs.machine.vatsettings.fastfield_cb;  % 0.1;
thresh = options.prefs.machine.vatsettings.fastfield_ethresh; % 0.2;

if useSI
    thresh=thresh.*(10^3);
end

if ~any(S.activecontacts{side}) % empty VAT, no active contacts.
    ofv.vertices=[0,0,0
        0,0,0
        0,0,0];
    ofv.faces=[1,2,3];
    varargout{1}=ofv;
    varargout{2}=0;
    return
end

resultfig=getappdata(lgfigure,'resultfig');
elstruct=getappdata(resultfig,'elstruct');
options=getappdata(resultfig,'options');
elspec=getappdata(resultfig,'elspec');
options.usediffusion=0;
coords=acoords{side};
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

options=ea_resolve_elspec(options);
Electrode_type = elspec.matfname;

Efield_all=zeros(100,100,100);
for source=S.sources
    stimsource=S.([sidec,'s',num2str(source)]);
    % constvol is 1 for constant voltage and 0 for constant current.
    amp1 = stimsource.amp;
    if amp1>0
        load([ea_getearoot,'templates',filesep,'standard_efields' filesep 'standard_efield_' Electrode_type '.mat']);
        count1=1;

        for cnt=1:length(cnts)
            % FastField only has monopolar (cathode) mode
            % So, if VC, all have 100%, no splitting
            if stimsource.(cnts{cnt}).pol==2
                ea_warndlg("Anodes are not supported in FastField")
                return
            end

            if stimsource.va==1
                if stimsource.(cnts{cnt}).perc ~= 0.0
                    perc(cnt) = 100;
                else
                    perc(cnt) = 0;
                end
            else
                perc(cnt) = stimsource.(cnts{cnt}).perc;
            end
            if perc(cnt)>0
                Im(count1)=stimsource.(cnts{cnt}).imp;
                count1=count1+1;
            end
        end

        constvol=stimsource.va==1;
        if constvol
            amp_mode = 'V';
            impedence = mean(Im)*1000;
        else
            amp_mode = 'mA';
            impedence = [];
        end

        [Efield2] = ea_get_efield(perc,standard_efield,amp1,conductivity,amp_mode,impedence);
        Efield_all = Efield_all+Efield2;
    end
end

Efield = Efield_all;

electrode_patient = elstruct;
load([ea_getearoot,'templates',filesep,'electrode_models',filesep,Electrode_type '.mat']);

[trans_mat,~,xg,yg,zg] = get_trans_mat(electrode,electrode_patient,grid_vec,side);

gv=grid_vec;

ea_dispt('Creating nifti header for export...');
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
Vvat.mat=mldivide([(chun1);(chun2);(chun3);ones(1,res)]',[gv{1}(chun1);gv{2}(chun2);gv{3}(chun3);ones(1,res)]')';
Vvat.mat = trans_mat * Vvat.mat;
Vvat.dim=[res,res,res];
Vvat.dt = [4, endian];
Vvat.n=[1 1];
Vvat.descrip='lead dbs - vat';

eeg = Efield;
eg=eeg;
eg=eg>thresh;
% binary e-field - "vat"

neeg=eeg;
neeg(~eg)=nan;

neeg(neeg>0)=ea_normal(neeg(neeg>0),1,0,' ',0,1,'TRUE');%
% normalized e-field (zscored).
neeg(~isnan(neeg))=neeg(~isnan(neeg))-min(neeg(~isnan(neeg)));
neeg(~isnan(neeg))=neeg(~isnan(neeg))/sum(neeg(~isnan(neeg))); % 0-1 distributed.

% [xg,yg,zg] = meshgrid(gv{1},gv{2},gv{3});

% eg=smooth3(eg,'gaussian',[25 25 25]);
ea_dispt('Calculating volume...');

for t=1:3
    spacing(t) = grid_vec{t}(2)-grid_vec{t}(1);
end

vatvolume=sum(eg(:))*spacing(1)*spacing(2)*spacing(3); % returns volume of vat in mm^3
S.volume(side)=vatvolume;

ea_dispt('Writing files...');

stimDir = fullfile(options.subj.stimDir, ea_nt(options), stimname);
ea_mkdir(stimDir);
filePrefix = ['sub-', options.subj.subjId, '_sim-'];

switch side
    case 1
        Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-fastfield_hemi-R.nii'];
        Vvate=Vvat;
        Vvatne=Vvat;
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-fastfield_hemi-R.nii'];
        Vvatne.fname = [stimDir, filesep, filePrefix, 'efieldgauss_model-fastfield_hemi-R.nii'];
    case 2
        Vvat.fname = [stimDir, filesep, filePrefix, 'binary_model-fastfield_hemi-L.nii'];
        Vvate = Vvat;
        Vvatne = Vvat;
        Vvate.fname = [stimDir, filesep, filePrefix, 'efield_model-fastfield_hemi-L.nii'];
        Vvatne.fname = [stimDir, filesep, filePrefix, 'efieldgauss_model-fastfield_hemi-L.nii'];
end

ea_savestimulation(S,options);

Vvate.img=eeg;
Vvate.dt = [16, endian];
ea_write_nii(Vvate);

Vvatne.img=neeg;
ea_write_nii(Vvatne);

Vvat.img=eg;
ea_write_nii(Vvat);

ea_dispt('Calculating isosurface to display...');
vatfv=isosurface(xg,yg,zg,Vvat.img,0.75);

caps=isocaps(xg,yg,zg,Vvat.img,0.5);

vatfv.faces=[vatfv.faces;caps.faces+size(vatfv.vertices,1)];
vatfv.vertices=[vatfv.vertices;caps.vertices];

try
    vatfv=ea_smoothpatch(vatfv,1,35);
catch
    try
        cd([ea_getearoot,'ext_libs',filesep,'smoothpatch']);
        mex ea_smoothpatch_curvature_double.c -v
        mex ea_smoothpatch_inversedistance_double.c -v
        mex ea_vertex_neighbours_double.c -v
        vatfv=ea_smoothpatch(vatfv);
    catch
        warndlg('Patch could not be smoothed. Please supply a compatible Matlab compiler to smooth VTAs.');
    end
end

% new save by Till to save VAT and quiver in seperate .mat-file for quick
% visualization
switch side
    case 1
        vatfvname = [stimDir, filesep, filePrefix, 'binary_model-fastfield_hemi-R.mat'];
    case 2
        vatfvname = [stimDir, filesep, filePrefix, 'binary_model-fastfield_hemi-L.mat'];
end

save(vatfvname,'vatfv','vatvolume');

% new vta.nii save, filled and eroded/dilated by 3 voxels.
Vvat.img=imfill(Vvat.img,'holes');
SE = strel('sphere',3);
Vvat.img = imerode(Vvat.img,SE);
Vvat.img = imdilate(Vvat.img,SE);
ea_write_nii(Vvat);

varargout{1}=vatfv;
varargout{2}=vatvolume;
ea_dispt('');
