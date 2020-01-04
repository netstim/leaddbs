function ea_concat_fMRI(directory,sessions,output)

if ~exist('output','var')
    output=[directory,'rest.nii'];
end
for sess=1:length(sessions)
   nii=ea_load_untouch_nii([directory,sessions{sess}]);
   if ~exist('allnii','var')
       allnii=nii;
       sessvec=ones(size(nii.img,4),1);
   else
       allnii.img=cat(4,allnii.img,nii.img);
       sessvec=[sessvec;repmat(sess,size(nii.img,4),1)];
   end
   
    
end

ea_save_untouch_nii(allnii,output);
save([directory,'rest_sessvec.mat'],'sessvec');


