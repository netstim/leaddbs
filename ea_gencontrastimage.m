function imat=ea_gencontrastimage(imat,howto)


switch howto
     case 2
        %% use raw tra version
        imat=squeeze(imat(:,:,:,1));

    case 3
        %% use smoothed tra only
        imat=squeeze(imat(:,:,:,1));

        imat = smooth3(imat,'gaussian',3);

    case 4
        %% use smoothed mean of tra and cor
        imat=mean(imat,4);

        imat = smooth3(imat,'gaussian',3);



    case 6
        %% use cor only
        imat=squeeze(imat(:,:,:,end));

   case 7
       %% use cor only but smooth it before.
        imat=squeeze(imat(:,:,:,end));

        imat = smooth3(imat,'gaussian',3);

    case 8
        %% use average of cor and tra:
        imat=mean(imat,4);
    case 9
        %% use cor*tra
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,end)));
    case 10
        %% use quadratic cor only
        imat=squeeze(imat(:,:,:,2).^4);
    case 11
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,end))).^4;
    case 12
        imat=double(squeeze(imat(:,:,:,1).*imat(:,:,:,end)));
        imat=imat.^4;

        imat = smooth3(imat,'gaussian',3);


    otherwise
        ea_error('Please set Contrast accurately.');

end

ea_delete('tmp.nii');
