function out=ea_query_sim_connectomematrix(seedmap,connectomename,outputfname)



s=ea_load_nii(seedmap);
seed=s.img(:);
out=s;
out.fname=outputfname;

C=matfile([ea_space([],'connectomes'),'fMRI',filesep,connectomename,filesep,'AllX.mat'],'Writable',false);
dataset=load([ea_space([],'connectomes'),'fMRI',filesep,connectomename,filesep,'dataset_volsurf.mat']);

ix=dataset.vol.outidx;

chunk=10000;
ea_dispercent(0,'Iterating connectome');
for vox=1:chunk:length(ix)
    try
        R=corr(seed(ix),double(C.X(1:length(ix),vox:vox+chunk-1)),'rows','pairwise');
        out.img(ix(vox:vox+chunk-1))=R;
    catch % last run:
        R=corr(seed(ix),double(C.X(1:length(ix),vox:length(ix))),'rows','pairwise');
        out.img(ix(vox:length(ix)))=R;
    end
    ea_dispercent(vox/length(ix));
    ea_write_nii(out);
end

ea_dispercent(1,'end');
ea_write_nii(out);
