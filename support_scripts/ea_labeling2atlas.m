function ea_labeling2atlas(labelname)
% this function will convert a whole-brain parcellation atlas (stored in
% /templates/labeling) into a Lead-DBS atlas (stored in /atlases).

suppress_suff=0; % cut this many characters from each atlas parcel, e.g. to delete "-R" or "-L" suffixes.

earoot=ea_getearoot;
odir=[ea_space([],'atlases'),labelname];
if ~exist(odir,'file')
    mkdir(odir)

    mkdir([odir,filesep,'lh'])
    mkdir([odir,filesep,'rh'])
    mkdir([odir,filesep,'midline'])
    mkdir([odir,filesep,'mixed'])
else
    ea_error('Atlas with this name already exists.');
end

fid=fopen([ea_space([],'labeling'),labelname,'.txt']);
A=textscan(fid,'%d %s');

parc=ea_load_nii([ea_space([],'labeling'),labelname,'.nii']);
parc.img=round(parc.img);

ea_dispercent(0,'Iterating atlas components');
cnt=1;

for p=A{1}'
    ea_dispercent(cnt/length(A{1}));

    thisp=parc.img;
    thisp(~(thisp==p))=0;
    thisp(thisp==p)=1;

    % determine whether left, right, mixed of midline

    [xx,yy,zz]=ind2sub(size(thisp),find(thisp));
    XYZ=[xx,yy,zz,ones(length(xx),1)]';
    XYZ=parc.mat*XYZ;
    hasleft=any(XYZ(1,:)<0);
    hasright=any(XYZ(1,:)>0);

    if hasleft && ~hasright % left only
        thisodir=[odir,filesep,'lh'];
    elseif ~hasleft && hasright % left only
        thisodir=[odir,filesep,'rh'];
    elseif hasleft && hasright % both hemispheres, either mixed or midline
        stats=bwconncomp(thisp);
        if stats.NumObjects>1 % mixed
            thisodir=[odir,filesep,'mixed'];
        elseif stats.NumObjects==1 % midline
            thisodir=[odir,filesep,'midline'];
        end
    end

    thisparcnii=parc;
    thisparcnii.img=thisp;
    thisparcnii.fname=[thisodir,filesep,A{2}{cnt}(1:end-suppress_suff),'.nii'];
        cnt=cnt+1;
    ea_write_nii(thisparcnii);
end
ea_dispercent(1,'end');
