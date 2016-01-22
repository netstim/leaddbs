atl=ea_load_nii('labels_on_colin_Nov2010.nii');
nucl=atl;
atl.img=round(atl.img);
labs={'LGN','MGN','Anterior nuclei','Central nuclei','Lateral dorsal nuclei','Lateral posterior nuclei','Medial dorsal nuclei','Pulvinar','VA','VL','VP'};
for nu=1:max(atl.img(:));
   nucl.img(:)=0;
   nucl.img(atl.img==nu)=1;
   if any(nucl.img(:))
   nucl.fname=[num2str(nu),'.nii'];
   spm_write_vol(nucl,nucl.img);
   end
end