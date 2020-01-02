function ea_roi2txt(roiname,txtoutput,upsamplefactor)
% used for DSI studio (exports .txt files identically as DSI studio does).
if ~exist('upsamplefactor','var')
    upsamplefactor=1;
end
if ~ismember(upsamplefactor,[1,2,4])
    ea_error('Only factors of 1, 2 or 4 are allowed.');
end
nii=ea_load_nii(roiname);
[xx,yy,zz]=ind2sub(size(nii.img),find(nii.img(:)));
xx=xx-1; % zero based vs. one based indexing
yy=yy-1; % zero based vs. one based indexing
zz=zz-1; % zero based vs. one based indexing
XYZ=[xx,yy,zz].*upsamplefactor;
n=size(XYZ,1);

oneXYZ=XYZ;
for x=0:(upsamplefactor-1)
    for y=0:(upsamplefactor-1)
        for z=0:(upsamplefactor-1)
            if any([x,y,z]) % exclude 0 0 0 case
                XYZ=[XYZ;...
                    oneXYZ+repmat([x y z],n,1)];
            end
        end
    end
end

fid=fopen(txtoutput,'w');
for entry=1:size(XYZ,1)
fprintf(fid,'%d %d %d\n',XYZ(entry,1),XYZ(entry,2),XYZ(entry,3));

end
fclose(fid);

