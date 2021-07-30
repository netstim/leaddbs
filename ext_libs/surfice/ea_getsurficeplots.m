function [im]=ea_getsurficeplots(niftiname,threshs)
% perform in tempdir given path handling issues in surfice
[pth,fn,ext]=fileparts(niftiname);
if ~exist(fullfile(pth,[fn,'_l_lat.png']),'file') || ~exist(fullfile(pth,[fn,'_l_med.png']),'file')
    tempdir=ea_getleadtempdir;
    uuid=ea_generate_uuid;
    tniftiname=[tempdir,uuid,'.nii'];
    copyfile(niftiname,tniftiname);

    [tpth,tfn,ext]=fileparts(tniftiname);
    ea_surficeoverlay(tniftiname,threshs,2);
    pause(0.5);
    delete(tniftiname);
    movefile(fullfile(tpth,[tfn,'_l_lat.png']),fullfile(pth,[fn,'_l_lat.png']));
    movefile(fullfile(tpth,[tfn,'_l_med.png']),fullfile(pth,[fn,'_l_med.png']));
end

imlat=imread(fullfile(pth,[fn,'_l_lat.png']));
immed=imread(fullfile(pth,[fn,'_l_med.png']));
im=fuse2im(imlat,immed);


function im=fuse2im(imlat,immed)

im=zeros(size(immed,1),size(imlat,2)+size(immed,2),3,'uint8');
im(1:size(immed,1),1:size(immed,2),1:3)=immed;
im(1:size(imlat,1),size(immed,2)+1:end,1:3)=imlat;
