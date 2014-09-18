function make_stn_czi_atlas


V=spm_vol([fileparts(which('lead')),filesep,'templates',filesep,'mni_hires_bb.nii']);
lcZI=spm_read_vols(V);
lSTN=spm_read_vols(V);
lcZI(:)=0;
lSTN(:)=0;
rcZI=spm_read_vols(V);
rSTN=spm_read_vols(V);
rcZI(:)=0;
rSTN(:)=0;

coords=[12.7 -14.6 -7.2 1
    -12.7 -14.6 -7.2 1
    12.6 -18.2 -5.8 1
    -12.6 -18.2 -5.8 1]'; % rSTN, lSTN, rcZI, lcZI. cZI after Blomstedt.

coords_vx=[V.mat\coords]';
coords_vx=round(coords_vx(:,1:3));

rSTN(coords_vx(1,1),coords_vx(1,2),coords_vx(1,3))=1;
lSTN(coords_vx(2,1),coords_vx(2,2),coords_vx(2,3))=1;
rcZI(coords_vx(3,1),coords_vx(3,2),coords_vx(3,3))=1;
lcZI(coords_vx(4,1),coords_vx(4,2),coords_vx(4,3))=1;


mkdir([fileparts(which('lead')),filesep,'atlases',filesep,'stnczi']);
mkdir([fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'lh']);
mkdir([fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'rh']);

V.fname=[fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'rh',filesep,'STN.nii'];
spm_write_vol(V,rSTN);
V.fname=[fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'lh',filesep,'STN.nii'];
spm_write_vol(V,lSTN);
V.fname=[fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'rh',filesep,'cZI.nii'];
spm_write_vol(V,rcZI);
V.fname=[fileparts(which('lead')),filesep,'atlases',filesep,'stnczi',filesep,'lh',filesep,'cZI.nii'];
spm_write_vol(V,lcZI);


