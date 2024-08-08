function ea_get_E_field_along_fibers(pt_folder,stim_space,e_field_file, MNI_connectome_file, side_suffix, threshold)
% get E-field metrics along fibers (computed in native space)
% By J.Roediger and K.Butenko

arguments
    pt_folder               % full path to the patient folder in the Lead-DBS dataset
    stim_space              % identifier for the stimulation ID and space, e.g. 'native/gs_111'
    e_field_file            % full path to the 4-D nifti (see ea_get_4Dfield_from_csv)
    MNI_connectome_file     % full path to the connectome, use merged_pathways.mat for pathway atlases
    side_suffix             % _rh - right hemisphere
    threshold               % sum and mean E-field metrics are only computed on segments above this E-field threshold (in V/m) 
end

% interpolation from the E-field grid to the fiber compartments
interpolation_method = 'linear';

% check which space is used
[space, stim_folder_name,~] = fileparts(stim_space);

% check the folder name if stored in data or merged_pathways.mat
[connectomePath,connectomeFileName,connectomeExtension] = fileparts(MNI_connectome_file);
if strcmp('data',connectomeFileName) || strcmp('merged_pathways',connectomeFileName)
    [~,connectomeName,~] = fileparts(connectomePath);
else
    ea_warndlg("MultiTract not supported, merge it. See ea_discfibers_merge_pathways.m")
    connectomeName = connectomeFileName;
end

pt_folder = char(pt_folder);
% warp the connectome to native if necessary
if strcmp(space,'native')
    connectomeFile_in_native = strcat([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, connectomeFileName, connectomeExtension]);
    % check if the connectome was already warped in native
    if isfile(connectomeFile_in_native)
        disp("The connectome was already warped, loading...")
        ftr = load(connectomeFile_in_native);
    else
        % otherwise warp the connectome to native space (and store in patient folder / connectomes)
        transform = ea_regexpdir([pt_folder, '/normalization/transformations'],'.*from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz',0,'f',0);
        anchor_img = ea_regexpdir([pt_folder,'/coregistration/anat'], '.*ses-preop_space-anchorNative.*',0,'f',0);
        ftr = ea_warp_fibers_MNI2native(pt_folder, MNI_connectome_file, transform{1}, anchor_img{1});
    end
else
    ftr = load(MNI_connectome_file);
end

% load the 4-D E-field nii 
% we assume that the first 3 components are X,Y,Z (see ea_get_4Dfield_from_csv)
nii = ea_load_nii(e_field_file);

%% Extract data from nii
efmat = nii.mat;
efdims = nii.dim;
for d = 1:3
    efield{d} = nii.img(:,:,:,d);           % Extract nii in x,y and z dimension
end

%% Get fibervectors and corresponding coordinates (mid), exclude points outside efield window

%disp("Computing E-Field projection along the fibers")

% exclude fiber coordinates which are outside of efield sample
bounds(:,1) = efmat*[1;1;1;1];
bounds(:,2) = efmat*[efdims';1];
bounds = bounds';
bounds(:,4) = [];

fibcoords = ftr.fibers(:,1:3);
fibidx = ftr.fibers(:,4);

outofbounds = any(fibcoords<=repmat(bounds(1,:),size(fibcoords,1),1,size(fibcoords,3))...
    |fibcoords>=repmat(bounds(2,:),size(fibcoords,1),1,size(fibcoords,3)),2);

fibcoords(outofbounds,:) = [];
fibidx(outofbounds,:) = [];

[allfibs,ci,ai] = unique(fibidx);

% exclude fibers with less than 3 coordinates
onefib = find(diff([ci;length(fibidx)])<3);
allfibs(onefib) = [];
[~,extmp] = ismember(onefib,ai);
fibidx(extmp) = [];
fibcoords(extmp,:) = [];

% Get fibervectors and corresponding coordinates
% check if there are any zeros
fibvector = diff(fibcoords);
fibmidcoords = fibcoords(1:end-1,:)+fibvector/2;
extmp = diff(fibidx)>0;
fibmidcoords(extmp,:) = [];
fibvector(extmp,:) = [];
fibmididx = fibidx(~[extmp;true]);

% Exclude fibers with large vector steps between the adjacent coordinates
thresh = 1.5;       % seems to be a good threshold for HCP atlas (template atlas has maximum values of 1)
extmp = any(abs(fibvector)>thresh,2);
% experc = sum(extmp)*100/length(fibmididx);
% 
% if experc<=1
%     disp([num2str(experc,2),' % of data points excluded'])
% else
%     warning([num2str(experc),...
%         ' % of data points excluded. This could point towards normalization issues'])
%     keyboard
% end

fibvector(extmp,:) = [];
fibmidcoords(extmp,:) = [];
fibmididx(extmp) = [];

fibmidcoords(:,4) = 1;

% Calculate efield vectors at fiber coordinates
efsub = fibmidcoords/efmat';
efsub(:,4) = [];
fibmidcoords(:,4) = [];

effiber = nan(size(fibvector));

if strcmp(interpolation_method,'neighbor')
    efsub = round(efsub);
    efind = sub2ind(efdims,efsub(:,1),efsub(:,2),efsub(:,3));
    
    for d = 1:3
        effiber(:,d) = efield{d}(efind);
    end
else
    % linear interpolation
     efsub_low = floor(efsub);
     efsub_high = ceil(efsub);
     efind_low = sub2ind(efdims,efsub_low(:,1),efsub_low(:,2),efsub_low(:,3));
     efind_high = sub2ind(efdims,efsub_high(:,1),efsub_high(:,2),efsub_high(:,3));
     relfc = mod(efsub,1);
     
     for d = 1:3
         eftmp_low = efield{d}(efind_low);
         eftmp_high = efield{d}(efind_high);
         effiber(:,d) = eftmp_low.*(1-relfc(:,d))+eftmp_high.*relfc(:,d);
     end
end

%% Obtain fiber E-metrics
% there could be nans in fibvector, because some connectomes have
% repetitions
fibvector_norm = fibvector./vecnorm(fibvector,2,2);
E_along_fiber = abs(dot(effiber,fibvector_norm,2));  % field projection onto the tangential unit vector, the directionality (the sign) is not relevant
E_magn_fiber = sqrt(dot(effiber,effiber,2));         % a magnitude of the E-field vector

% there is no sense to restore dimensionality of ftr, just keep reference to orig idx
fibers_E_proj = zeros(size(fibmididx,1),6);
fibers_E_proj(:,1:3) = fibmidcoords;  % we should recompute them in MNI?
fibers_E_proj(:,4) = fibmididx;
fibers_E_proj(:,5) = E_along_fiber;
fibers_E_proj(:,6) = E_magn_fiber;   % just for consistency

% we should consider different thresholds for these
% E-field projections on fiber
E_metrics.proj_peak = zeros(size(ftr.idx));
E_metrics.proj_5perc_peak = zeros(size(ftr.idx));
E_metrics.proj_sum = zeros(size(ftr.idx));   % only where E-metric > threshold
E_metrics.proj_mean = zeros(size(ftr.idx));  % only where E-metric > threshold

% E-field magnitudes on fiber
E_metrics.magn_peak = zeros(size(ftr.idx));
E_metrics.magn_5perc_peak = zeros(size(ftr.idx));
E_metrics.magn_sum = zeros(size(ftr.idx));   % only where E-metric > threshold
E_metrics.magn_mean = zeros(size(ftr.idx));  % only where E-metric > threshold

% we do not need to do any VAT intersections anymore, because we have all
% metrics available on fibers already

for fib_idx = 1:size(ftr.idx,1)
    fib_E_proj = fibers_E_proj(fibers_E_proj(:,4) == fib_idx,:);
    if ~isempty(fib_E_proj)

        E_metrics.proj_peak(fib_idx) = ea_nanmax(fib_E_proj(:,5));
        E_metrics.proj_5perc_peak(fib_idx) = ea_nanmean(maxk(fib_E_proj(:,5),ceil(0.05*numel(size(ftr.idx,1)))));

        % only compute sum and mean for segments above the threshold
        idx_thresh_proj =  find(fib_E_proj(:,5)>=threshold*0.001);
        if ~isempty(idx_thresh_proj)
            E_metrics.proj_sum(fib_idx) = ea_nansum(fib_E_proj(idx_thresh_proj,5));
            E_metrics.proj_mean(fib_idx) = ea_nanmean(fib_E_proj(idx_thresh_proj,5));
        end
    
        E_metrics.magn_peak(fib_idx) = ea_nanmax(fib_E_proj(:,6));
        E_metrics.magn_5perc_peak(fib_idx) = ea_nanmean(maxk(fib_E_proj(:,6),ceil(0.05*numel(size(ftr.idx,1)))));

        % only compute sum and mean for segments above the threshold
        idx_thresh_magn = find(fib_E_proj(:,6)>=threshold*0.001);
        if ~isempty(idx_thresh_magn)
            E_metrics.magn_sum(fib_idx) = ea_nansum(fib_E_proj(idx_thresh_magn,6));
            E_metrics.magn_mean(fib_idx) = ea_nanmean(fib_E_proj(idx_thresh_magn,6));
        end
    end
end

%% Store the projections in patient_folder/connectomes/dMRI/connectome_name/stim_folder_side
stim_folder_name = [stim_folder_name, side_suffix];

if ~isfolder([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, stim_folder_name])
    mkdir([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, stim_folder_name])
end

if strcmp(space,'native')
    save([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, stim_folder_name, filesep, 'E_metrics.mat'], 'E_metrics')
else
    ea_mkdir([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, stim_folder_name, filesep, 'MNI'])
    save([pt_folder, filesep,'connectomes',filesep,'dMRI', filesep, connectomeName, filesep, stim_folder_name, filesep,'MNI',filesep,'E_metrics.mat'], 'E_metrics')
end
