function kb2_get_E_field_along_fibers(pt_folder,stim_folder,e_field_file, MNI_connectome_file, side_suffix)

% get E-field along fibers (computed in native space)
% By J.Roediger and K.Butenko

% first warp the connectome to native space (and store in patient folder / miscellaneous)
[connectomePath,connectomeFileName,connectomeExtension] = fileparts(MNI_connectome_file);

% check the folder name if stored in data or merged_pathways.mat
if strcmp('data',connectomeFileName) || strcmp('merged_pathways',connectomeFileName)
    [~,connectomeName,~] = fileparts(connectomePath);
else
    ea_warndlg("MultiTract not supported, merge it. See ea_discfibers_merge_pathways.m")
    connectomeName = connectomeFileName;
end

pt_folder = char(pt_folder);

connectomeFile_in_native = strcat([pt_folder, filesep,'miscellaneous', filesep, connectomeName, filesep, connectomeFileName, connectomeExtension]);
% check if the connectome was already warped in native
if isfile(connectomeFile_in_native)
    disp("The connectome was already warped, loading...")
    ftr = load(connectomeFile_in_native);
else
    % modify to search for non-ants as well
    transform = ea_regexpdir([pt_folder, '/normalization/transformations'],'.*from-anchorNative_to-MNI152NLin2009bAsym_desc-ants.nii.gz',0,'f',0);
    anchor_img = ea_regexpdir([pt_folder,'/coregistration/anat'], '.*ses-preop_space-anchorNative.*',0,'f',0);
    ftr = ea_warp_fibers_MNI2native(pt_folder, MNI_connectome_file, transform{1}, anchor_img{1});
end
% load the 4-D E-field nii (in native space!)
% , we assume that the first 3 components are X,Y,Z
nii = ea_load_nii(e_field_file);

% interpolation_method = 'neighbor';
interpolation_method = 'linear';
    
%% Extract data from nii
efmat = nii.mat;
efdims = nii.dim;
for d = 1:3
    %efield{d} = nii.img(:,:,:,d+1);        % Extract nii in x,y and z dimension
    efield{d} = nii.img(:,:,:,d);        % Extract nii in x,y and z dimension
end

%% Get fibervectors and corresponding coordinates (mid), exclude points outside efield window

disp("Computing E-Field projection along the fibers")

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


%% Obtain fiber orientation in relation to efield vector
fibvector_norm = fibvector./vecnorm(fibvector,2,2);
E_along_fiber = abs(dot(effiber,fibvector_norm,2));

% there is no sense to restore dimensionality of ftr, just keep reference to orig idx
fibers_E_proj = zeros(size(fibmididx,1),5);
fibers_E_proj(:,1:3) = fibmidcoords;  % we should recompute them in MNI?
fibers_E_proj(:,4) = fibmididx;
fibers_E_proj(:,5) = E_along_fiber;

% we can also compute 5% peak here, etc.
% But better pass the whole thing to FF
E_peak = zeros(size(ftr.idx));
E_5perc_peak = zeros(size(ftr.idx));

for fib_idx = 1:size(ftr.idx,1)
    fib_E_proj = fibers_E_proj(fibers_E_proj(:,4) == fib_idx,:);
    if isempty(fib_E_proj)
        E_peak(fib_idx) = 0.0;
        E_5perc_peak(fib_idx) = 0.0;
    else
        E_peak(fib_idx) = max(fib_E_proj(:,5));
        E_5perc_peak(fib_idx) = mean(maxk(fib_E_proj(:,5),5));
    end
end
    
%% Store the projections in patient_folder/miscellaneous/connectome_name/stim_folder_side
[~,stim_folder_name,~] = fileparts(stim_folder);
stim_folder_name = [stim_folder_name, side_suffix];

if ~isfolder([pt_folder, filesep,'miscellaneous', filesep, connectomeName, filesep, stim_folder_name])
    mkdir([pt_folder, filesep,'miscellaneous', filesep, connectomeName, filesep, stim_folder_name])
end

save([pt_folder, filesep,'miscellaneous', filesep, connectomeName, filesep, stim_folder_name, filesep, 'E_peak.mat'], 'E_peak')
save([pt_folder, filesep,'miscellaneous', filesep, connectomeName, filesep, stim_folder_name, filesep, 'E_5perc_peak.mat'], 'E_5perc_peak')
