function ea_ftr2nii(ftr,niftifile)


if ~exist('niftifile','var')
    niftifile=[ea_space,'t1.nii'];
end
nii=ea_load_nii(niftifile);

[fibers,idx,voxmm,mat,vals]=ea_loadfibertracts(ftr);

if strcmp(voxmm,'vox') % convert to mm based on mat supplied with ftr file
    T=(mat*[fibers(:,1:3)';ones(1,size(fibers,1))]);
    fibers(:,1:3)=round(T(1:3,:)');
end

T=(nii.mat\[fibers(:,1:3)';ones(1,size(fibers,1))]);
fibers=round(T(1:3,:)');

nii.img(:)=0;
niisz=size(nii.img);

fibers(fibers<=0)=1;
for dim=1:3
fibers(fibers(:,dim)>niisz(dim),dim)=niisz(dim);
end
cnt=1;
ea_dispercent(0,'Iterating fibers');
for fib=1:length(idx)
    if ~isnan(vals(fib)) && ~isinf(vals(fib))
        fixx=cnt:cnt+(idx(fib)-1);
        ixx=sub2ind(niisz,fibers(fixx,1),fibers(fixx,2),fibers(fixx,3));
        nii.img(ixx)=nii.img(ixx)+vals(fib);
    end
    cnt=cnt+idx(fib);
    ea_dispercent(fib/length(idx));
end
ea_dispercent(1,'end');

[pth,fn,ext]=fileparts(ftr);
nii.fname=fullfile(pth,[fn,'.nii']);
nii.dt=[16,0];
ea_write_nii(nii);
