function [dtd err] = nifti_to_DTD(data_names,T)

tensor = nifti_to_mrstruct('series3D',data_names);
evals = zeros(size(tensor.dataAy,1),size(tensor.dataAy,2),size(tensor.dataAy,3),3);
evecs = zeros(size(tensor.dataAy,1),size(tensor.dataAy,2),size(tensor.dataAy,3),3,3);

for x = 1:size(tensor.dataAy,1),
    for y = 1:size(tensor.dataAy,2),
        for z = 1:size(tensor.dataAy,3),
            D = [tensor.dataAy(x,y,z,1) tensor.dataAy(x,y,z,2) tensor.dataAy(x,y,z,3);
                 tensor.dataAy(x,y,z,2) tensor.dataAy(x,y,z,4) tensor.dataAy(x,y,z,5);
                 tensor.dataAy(x,y,z,3) tensor.dataAy(x,y,z,5) tensor.dataAy(x,y,z,6)];
             if exist('T'),
                 D = T*D*T';
             end;
             if any(not(isnan(D(:)))),
                 [U d] = eig(D);
                 evals(x,y,z,:) = diag(d);
                 evecs(x,y,z,:,:) = U;
             end;
        end;
    end;
end;

evecs(:,:,:,[1 2 3],:) = evecs(:,:,:,[2 1 3],:);
Eval = mrstruct_init('series3D',evals,tensor);
Evec = mrstruct_init('series3DEchos',evecs,tensor);
if size(tensor.dataAy) > 6,
    b0 = mrstruct_init('volume',mean(tensor.dataAy(:,:,:,7:end),4),tensor);    
else,
    b0 = mrstruct_init('volume',evals(:,:,:,1)*0+1,tensor);    
end;

[dtd err] = dtdstruct_init('DTD',Evec,Eval,b0,'b0_image_struc');

return;

