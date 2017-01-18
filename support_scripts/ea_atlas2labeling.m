function ea_atlas2labeling(atlasname,spacefile)


copyfile([ea_space([],'atlases'),atlasname],[ea_space([],'atlases'),atlasname,'_copy']);
if ~exist('spacefile','var')
spacefile=[ea_space,'t2.nii'];
end

copybase=[ea_space([],'atlases'),atlasname,'_copy',filesep];

subfs={'lh','rh','mixed','midline'};
cnt=1;
for subf=1:length(subfs)
    subbase=[copybase,subfs{subf},filesep];
    switch subfs{subf}
        case 'lh'
            append='_L';
        case 'rh'
            append='_R';
        otherwise
            append='';
    end

   system(['gunzip ',subbase,'*.nii.gz']);

   nuclei=dir([subbase,'*.nii']);

   for nucleus=1:length(nuclei)
       ea_conformspaceto(spacefile,[subbase,nuclei(nucleus).name],1);
       [~,nuclname]=fileparts(nuclei(nucleus).name);
       Astr{cnt}=[nuclname,append];
       nii=ea_load_nii([subbase,nuclei(nucleus).name]);
       if ~exist('Avol','var') % init
           Avol=zeros(size(nii.img));
       end
       Avol(:,:,:,cnt)=nii.img;
       cnt=cnt+1;
   end

end

Avol(isnan(Avol))=0;
Avol(Avol<0.5)=0;
[~,Fvol]=max(Avol,[],4);
lab=nii;
lab.img=Fvol;
lab.fname=[ea_space([],'labeling'),atlasname,'.nii'];
ea_write_nii(lab);

f=fopen([ea_space([],'labeling'),atlasname,'.txt'],'w');
for nucl=1:length(Astr)
    fprintf(f,'%d %s\n',nucl,Astr{nucl});
end
fclose(f);

rmdir([ea_space([],'atlases'),atlasname,'_copy'],'s')











