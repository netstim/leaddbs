function ringRemoval(niinames)

for k = 1:length(niinames),
    
    [pa na ext] = fileparts(niinames{k});
    if strcmp(ext,'.nii'),    
        im = nifti(niinames{k});
        for i = 1:size(im.dat,4),
            v = double(im.dat(:,:,:,i));
            v = ringRm(v,[1 3 10]);
            im.dat(:,:,:,i) = v;
            fprintf('.');
        end;
    elseif strcmp(ext,'.mat');
        mr = mrstruct_read(niinames{k});
        for i = 1:size(mr.dataAy,4),
            v = double(mr.dataAy(:,:,:,i));
            v = ringRm(v,[1 3 10]);
            mr.dataAy(:,:,:,i) = v;
            fprintf('.');
        end;
        mrstruct_write(mr,niinames{k});
    else
        warning('filetype not recognized');
    end
        
end;
fprintf('\n');
