%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  dtd = convertHARDI2DTD(hr)
%
%  converts a HARDI mrstruct into a DTD
%
%  hr - a hardi mrstruct
%  dtd - the resulting tensor file
%
%  in order to save the dtd >> dtstruct_write(dtd,'mynew_DTD.mat')
%
%

function dtd = ea_convertHARDI2DTD(hr,threshold,kurto)

if not(exist('kurto')),
    kurto = true;
end;


if not(exist('threshold')),
    threshold = 0;
end;

if isempty('threshold'),
    threshold = 0;
    
end;

signal = hr.dataAy;
bTensor = hr.user.bTensor;
vox = hr.vox;
edges = hr.edges;



sz = size(signal);



mean_DTI = mean(signal(:,:,:,:),4);

display('calculating diffusion tensor ...');

Eval = zeros([sz(1:3) 3]);
Evec = zeros([sz(1:3) 3 3]);
Err = zeros(sz(1:3));
Kval = zeros([sz(1:3) 8]);


bval = squeeze(bTensor(1,1,:)+bTensor(2,2,:)+bTensor(3,3,:));
b0idx = find(bval<100);
b0_image = mean(signal(:,:,:,b0idx),4);


for s = 1:size(signal,3),
    [eval evec err sqerr kurtParams] = ea_DWI_calculation(bTensor, squeeze(signal(:,:,s,:)), [sz(1), sz(2)],threshold,kurto);
    Eval(:,:,s,:) = eval;
    Evec(:,:,s,:,:) = evec;
    Err(:,:,s) = err;   
    if not(isempty(kurtParams)),
        Kval(:,:,s,:) = kurtParams;        
    else
        kurto = false;
    end;
end;


display('saving ...');


%% save as dtd_struct
b0_image_struc= hr;
b0_image_struc.dataAy = b0_image;
b0_image_struc.dim4 = 'unused';
b0_image_struc.memoryType = 'volume';



meanDWI_image_struc = mrstruct_init('volume',mean_DTI,b0_image_struc);
EigenVect = mrstruct_init('series3DEchos',Evec,b0_image_struc);
EigenVal = mrstruct_init('series3D',Eval,b0_image_struc);
if kurto,
    
    K_axial = mrstruct_init('volume',Kval(:,:,:,1),b0_image_struc);
    K_orth = mrstruct_init('volume',Kval(:,:,:,2),b0_image_struc);
    K_mean = mrstruct_init('volume',Kval(:,:,:,3),b0_image_struc);
    D_orth = mrstruct_init('volume',Kval(:,:,:,4),b0_image_struc);
    D_para_int = mrstruct_init('volume',Kval(:,:,:,5),b0_image_struc);
    D_para_ext = mrstruct_init('volume',Kval(:,:,:,6),b0_image_struc);
    D_orth_ext = mrstruct_init('volume',Kval(:,:,:,7),b0_image_struc);
    volfrac = mrstruct_init('volume',Kval(:,:,:,8),b0_image_struc);
    
    [dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, b0_image_struc, 'b0_image_struc',meanDWI_image_struc,'meanDWI_image_struc',...
                                        K_axial,'K_axial', K_orth ,'K_orth_maxK', K_mean, 'K_mean',D_para_int ,'D_para_int' ,D_para_ext,'D_para_ext', ...
                                         D_orth_ext , 'D_orth_ext' , volfrac,'volfrac_int');
                                    
else
    [dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, b0_image_struc, 'b0_image_struc',meanDWI_image_struc,'meanDWI_image_struc');
end;
                                    
                                    
if ~isempty(errStr)
    error(errStr)
end



% sort eigenvals and eigenvects
[dtd, errStr] = dtdstruct_modify(dtd,'sortEigvec');
if ~isempty(errStr)
    error(errStr)
end


return;








