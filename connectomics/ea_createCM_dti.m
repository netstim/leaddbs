function [DTI_CM, DTI_LEN] = ea_createCM_dti(options)
% This function creates a structural (dMRI-based) Connectivity matrix.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% Use only fiber endpoints for structural CM.
% This increases the Chance, dass Seeds/Terminals in (grauer) GM mit
% Atlas-Labels fallen, insbesondere bei fokalen DWI-FOVs.
useendpointsonly = 1;
ea_warp_parcellation(options.prefs.b0, options);

minlen   = options.prefs.lc.struc.minlen;
directory = [options.root, options.patientname, filesep];

%% get node definition of current parcellation scheme
Vatl = ea_load_nii([directory, 'templates', filesep, 'labeling', filesep, ...
                    'b0w', options.lc.general.parcellation, '.nii,1']);

%% get fiber definition
disp('Loading FTR-File.');
[fibs, idx] = ea_loadfibertracts([options.root, options.patientname, filesep, options.prefs.FTR_unnormalized]);
fibs = double(fibs);

%% create CM
disp('Initializing structural CM.');

aID = fopen([ea_space(options, 'labeling'), options.lc.general.parcellation, '.txt']);
atlas_lgnd = textscan(aID, '%d %s');
d = length(atlas_lgnd{1}); % number of ROIs
DTI_CM  = zeros(d);
DTI_LEN = zeros(d);

if useendpointsonly
    [DTI_CM, DTI_LEN] = ea_createCM_dti_endpoints(fibs, idx, Vatl, atlas_lgnd{1}, minlen, DTI_CM, DTI_LEN);
else
    DTI_CM = ea_createCM_dti_tracts(fibs, idx, Vatl, atlas_lgnd{1}, minlen, DTI_CM);
end


function DTI_CM = ea_createCM_dti_tracts(fibs, idx, Vatl, atlasIndices, minlen, DTI_CM)

% Sample atlas labels along fibers using the current coordinate convention.
% In the classic Lead-Connectome pipeline, fibers are assumed to be in
% mm/world space compatible with Vatl.mat. However, depending on how
% FTR.mat was generated (legacy vs. BIDS, voxel vs. mm storage), this
% assumption may be violated, leading to effectively zero atlas hits and
% thus an empty connectivity matrix. To make the pipeline robust across
% both conventions, we (1) try the classic mm-based sampling and (2) if we
% observe suspiciously few non-zero atlas labels, we fall back to treating
% fiber coordinates as voxel indices and explicitly mapping them to mm
% space via Vatl.mat before sampling again. Whichever variant yields more
% atlas hits is used.

% Primary: assume fibs are in mm/world space
fib2parc = round(spm_sample_vol(Vatl, fibs(:,1), fibs(:,2), fibs(:,3), 0));

% Diagnostics: check how many points fall into non-zero atlas labels
nonzero_primary = sum(fib2parc(:) ~= 0 & ~isnan(fib2parc(:)));
total_samples   = numel(fib2parc);

% If we have almost no atlas hits, try voxel->mm fallback
hit_threshold = max(10, round(0.001 * total_samples)); % at least 10 or 0.1%
if nonzero_primary < hit_threshold
    fprintf(['WARNING ea_createCM_dti: Very few atlas hits in primary sampling ', ...
             '(%d/%d non-zero). Trying voxel->mm fallback mapping...\n'], ...
            nonzero_primary, total_samples);

    pts_vox = [fibs(:,1:3), ones(size(fibs,1),1)]';
    pts_mm  = Vatl.mat * pts_vox;
    fib2parc_fb = round(spm_sample_vol(Vatl, pts_mm(1,:)', pts_mm(2,:)', pts_mm(3,:)', 0));
    nonzero_fb = sum(fib2parc_fb(:) ~= 0 & ~isnan(fib2parc_fb(:)));

    if nonzero_fb > nonzero_primary
        fprintf(['INFO ea_createCM_dti: Using voxel->mm fallback mapping ', ...
                 '(%d/%d non-zero samples).\n'], nonzero_fb, total_samples);
        fib2parc = fib2parc_fb;
    else
        fprintf(['INFO ea_createCM_dti: Fallback mapping did not improve atlas hits ', ...
                 '(primary %d, fallback %d). Keeping primary mapping.\n'], ...
                nonzero_primary, nonzero_fb);
    end
end
cnt       = 1;
fibercount = length(idx);
ea_dispercent(0, ['Iterating through ', num2str(fibercount), ' fibers']);

for fiber = 1:fibercount
    if idx(fiber) > minlen % only include fibers > minimum length
        % find all the regions the fiber goes through
        thisfibconnects = unique(fib2parc(cnt:cnt+idx(fiber)-1));
        % map atlas indices (the indices in the atlas may not be from 1 to atlas size)
        thisfibconnects = find(ismember(atlasIndices, thisfibconnects));
        thisfibconnects = thisfibconnects(thisfibconnects > 0);
        % locate the regions in the connectivity matrix (diagonal values are also set)
        conmesh    = meshgrid(thisfibconnects, thisfibconnects);
        matindices = sub2ind(size(DTI_CM), conmesh, conmesh');

        DTI_CM(matindices) = DTI_CM(matindices) + 1;
        cnt = cnt + idx(fiber);
        ea_dispercent(fiber / fibercount);
    end
end

ea_dispercent(1, 'end');


function [DTI_CM, DTI_LEN] = ea_createCM_dti_endpoints(fibs, idx, Vatl, atlasIndices, minlen, DTI_CM, DTI_LEN)

disp('Calculating seeds and terminals...');
fibercount = length(idx);
seeds = zeros(fibercount, 3);
terms = zeros(fibercount, 3);

cnt = 1;
for fiber = 1:fibercount
    seeds(fiber,:) = fibs(cnt, 1:3);
    terms(fiber,:) = fibs(cnt+idx(fiber)-1, 1:3);
    cnt = cnt + idx(fiber);
end

if any(seeds(:) < -1) || any(terms(:) < -1) % mm-notation, convert to voxel notation
   seeds = [seeds(:,1), seeds(:,2), seeds(:,3), ones(size(seeds,1),1)]';
   terms = [terms(:,1), terms(:,2), terms(:,3), ones(size(terms,1),1)]';
   seeds = (Vatl.mat \ seeds)';
   terms = (Vatl.mat \ terms)';
   seeds = seeds(:,1:3);
   terms = terms(:,1:3);
end

seedIDX = round(spm_sample_vol(Vatl, seeds(:,1), seeds(:,2), seeds(:,3), 0));
termIDX = round(spm_sample_vol(Vatl, terms(:,1), terms(:,2), terms(:,3), 0));

vizz = 0;
if vizz
    figure;
    plot3(seeds(:,1), seeds(:,2), seeds(:,3), 'r.');
    nii = ea_load_nii(Vatl.fname);
    hold on
    [xx, yy, zz] = ind2sub(size(nii.img), find(nii.img(:) > 0));
    plot3(xx, yy, zz, 'g.');
end

clear seeds terms

ea_dispercent(0, ['Iterating through ', num2str(length(idx)), ' fibers']);
conns = 0;

endpts = [seedIDX, termIDX];
conIDX = endpts > 0;
conIDX = sum(conIDX, 2);
conIDX = find(conIDX == 2); % only take fibers into account that dont start/end in zero-regions
for fiber = conIDX'
    percent = fiber / fibercount;
    ea_dispercent(percent);
    if idx(fiber) > minlen % only include fibers >minimum length
        %if sum(thisfib(3,:)<-55)<1 % check if fiber exits the brain through spinal chord.. choose a cutoff(i.e.=-55mm).
        DTI_CM(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber)) =  ...
            DTI_CM(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber)) + 1;

        DTI_CM(atlasIndices == termIDX(fiber), atlasIndices == seedIDX(fiber)) =  ...
            DTI_CM(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber));  % symmetrize Matrix.

        DTI_LEN(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber)) =  ...
            DTI_LEN(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber)) + idx(fiber);

        DTI_LEN(atlasIndices == termIDX(fiber), atlasIndices == seedIDX(fiber)) =  ...
            DTI_LEN(atlasIndices == seedIDX(fiber), atlasIndices == termIDX(fiber));  % symmetrize Matrix.

        conns = conns + 1; % connection count
        %end
    end
end

DTI_LEN = DTI_LEN ./ DTI_CM;
DTI_LEN(isnan(DTI_LEN) | isinf(DTI_LEN)) = 0;

ea_dispercent(1, 'end');

disp(['In total used ', num2str(conns), '/', num2str(fiber), ...
      ' fibers to connect ', num2str(length(DTI_CM)), ' regions.']);

disp('Done.');

