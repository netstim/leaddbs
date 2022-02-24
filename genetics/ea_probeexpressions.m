function sims=ea_probeexpressions(testmap,db,generange,smoothkernel)

if exist('db','var') && ~isempty(db)
    load(db);
else
    prefs=ea_prefs;
    load([prefs.genetics.dbdir,'genedb.mat']);
end
if ~exist('smoothkernel','var')
    smoothkernel=0;
end
if ~exist('generange','var')
    generange=1:size(genedb,2);
end
if isempty(generange)
    generange=1:size(genedb,2);
end

if smoothkernel
    uid=ea_generate_uuid;
    copyfile(testmap,[ea_getleadtempdir,uid,'.nii']);
    spm_smooth([ea_getleadtempdir,uid,'.nii'],[ea_getleadtempdir,'s',uid,'.nii'],[smoothkernel smoothkernel smoothkernel]);
    testnii=ea_open_vol([ea_getleadtempdir,'s',uid,'.nii']);
else
    testnii=ea_open_vol(testmap);
end

if size(idxmm,1)==3
    idxmm=[idxmm;ones(1,size(idxmm,2))];
end

querypoints=testnii.mat\idxmm; % in voxel space of test map image (which has to be in MNI space)

vals=spm_sample_vol(testnii,querypoints(1,:),querypoints(2,:),querypoints(3,:),1)';

sims=corr(vals,single(genedb(:,generange)),'rows','pairwise');


% figure;
% hist(Rs,1000)


