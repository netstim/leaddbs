function imat=gencontrastimage(imat,howto)

switch howto
    case 1
        %% use cor only
        imat=squeeze(imat(:,:,:,2));

   case 2
       %% use cor only but smooth it before.
        imat=squeeze(imat(:,:,:,2));
        nii=make_nii(imat);
        save_nii(nii,'tmp.nii');
        spm_smooth('tmp.nii','tmp.nii',[3 3 3],0);
        nii=load_nii('tmp.nii');
        imat=nii.img;
    case 3
        %% use average of cor and tra:
        imat=mean(imat,4);
    case 4
        %% use cor*tra
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,2)));
    case 5
        %% use quadratic cor only
        imat=squeeze(imat(:,:,:,2).^4);
    case 6
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,2))).^4;
    case 7
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,2)));
        imat=imat.^4;
        nii=make_nii(imat);
        save_nii(nii,'tmp.nii');
        spm_smooth('tmp.nii','tmp.nii',[3 3 3],0);
        nii=load_nii('tmp.nii');
        imat=nii.img;
        
    case 8
        %% use smoothed tra only
        imat=squeeze(imat(:,:,:,1));
        nii=make_nii(imat);
        save_nii(nii,'tmp.nii');
        spm_smooth('tmp.nii','tmp.nii',[3 3 3],0);
        nii=load_nii('tmp.nii');
        imat=nii.img;
        
  case 9
        %% use smoothed mean of tra and cor
        imat=mean(imat,4);
        nii=make_nii(imat);
        save_nii(nii,'tmp.nii');
        spm_smooth('tmp.nii','tmp.nii',[3 3 3],0);
        nii=load_nii('tmp.nii');
        imat=nii.img;
        
    case 10
        %% use raw tra version
                imat=squeeze(imat(:,:,:,1));
        
end

if exist('tmp.nii','file')
delete('tmp.nii')
end