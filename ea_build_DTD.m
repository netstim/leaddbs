function [mr,b0idx,finames] = ea_build_DTD(bvalue,dirname,dtif,DTDf,HARDIf,bvalf,bvecf)
% converts the dti.nii and dti.bval/.bvec files to Freiburg DTD and raw
% HARDI-formats.


V=spm_vol([dirname,dtif]);

for fi=1:length(V)
   finames{fi}=[dirname,dtif,',',num2str(fi)]; 
end




%% change into mrstruct format

mr =  ea_nifti_to_mrstruct('series3D',finames);


sz = size(mr.dataAy);

%% build dtd

mr.user.bfactor = bvalue; % enter bvalue


[mr.user.bDir, b0idx] = grads(dirname,bvecf,bvalf); % gradient directions



for k = 1:size(mr.user.bDir,2),
    mr.user.bTensor(:,:,k) = mr.user.bDir(:,k)*mr.user.bDir(:,k)';
end;

b0_image = mean(mr.dataAy(:,:,:,b0idx),4);
mean_DTI = mean(mr.dataAy(:,:,:,setdiff(1:sz(4),b0idx)),4);

display('Concatenating HARDI-Data ...');

Eval = zeros([sz(1:3) 3]);
Evec = zeros([sz(1:3) 3 3]);
Err = zeros(sz(1:3));

for s = 1:size(mr.dataAy,3),
    [eval evec err] = do_calculation(mr.user.bTensor, mr.user.bfactor, squeeze(mr.dataAy(:,:,s,:)), [size(mr.dataAy,1),size(mr.dataAy,2)],[1:7],0);
    Eval(:,:,s,:) = eval;
    Evec(:,:,s,:,:) = evec;
    Err(:,:,s) = err;   
end;



display('Done. Saving ...');


%% save as dtd_struct
b0_image_struc= mrstruct_init('volume', b0_image);
b0_image_struc.tr = mr.tr;
b0_image_struc.te = mr.te;
b0_image_struc.ti = mr.ti;
b0_image_struc.patient = mr.patient;
b0_image_struc.orient = mr.orient;
b0_image_struc.edges = mr.edges;
b0_image_struc.user.dti_calc_Ver = 'version_xxx';
b0_image_struc.vox = mr.vox;

meanDWI_image_struc = mrstruct_init('volume',mean_DTI,b0_image_struc);
EigenVect = mrstruct_init('series3DEchos',Evec,b0_image_struc);
EigenVal = mrstruct_init('series3D',Eval,b0_image_struc);
[dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, b0_image_struc, 'b0_image_struc',meanDWI_image_struc,'meanDWI_image_struc');


dtdname = [dirname DTDf];
hardiname = [dirname HARDIf];


% sort eigenvals and eigenvects
dtd = dtdstruct_modify(dtd,'sortEigvec');

% write dtd
[res, errStr] = dtdstruct_write(dtd, dtdname);
if isempty(errStr)
    msgStr= sprintf('DTI calculations done, write file dtdStruct as %s',dtdname);
else
    ea_error('Error: DTI calculations done, BUT could not write DTD');
    msgStr= sprintf('Error: DTI calculations done, BUT could not write DTD');
end

% write HARDI
 
[res] = mrstruct_write(mr, hardiname);












return;




% function do_calculation, copied from ika, but shortened and cleaned
function [EigenVal, eigVect, M_error, sqErr] = do_calculation(DEscheme, bfactor, singleSliceDWI_set, matrix,indB0s,threshold)

%% calculation of the b-matrix
%check if already bmatrix (necessary for bruker)
DEscheme = squeeze(DEscheme);
DE_size = size(DEscheme);
if DE_size(1) == 3 && DE_size(2) == 3
    B_tensor = DEscheme;
else
    B_tensor =zeros(3,3,length(DEscheme));

    for i = 1:DE_size(1)
        tmp = DEscheme(i,:)' * DEscheme(i,:);
        if sum(sum(abs(tmp) )) > 0 % ika 041118
            B_tensor(:,:,i)= (bfactor/trace(tmp)) * tmp;
        end
    end
end
nsteps =length(B_tensor);
X1 = [reshape(B_tensor(1,1,:),1,nsteps,1)]; %1st diag element
X2 = [reshape(B_tensor(2,2,:),1,nsteps,1)]; %2nd diag element
X3 = [reshape(B_tensor(3,3,:),1,nsteps,1)]; %3rd diag element
X4 = [reshape(B_tensor(1,2,:),1,nsteps,1)];
X5 = [reshape(B_tensor(1,3,:),1,nsteps,1)];
X6 = [reshape(B_tensor(2,3,:),1,nsteps,1)];
anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
iAnisoX=pinv(anisoX);

%%   do dti calculation for the particular slice
b0Mean= mean(singleSliceDWI_set(:,:,indB0s), 3);
mask= double(b0Mean > threshold);


%% initiate variables
eigVect=zeros(matrix(1), matrix(2), 3,3);
difComp=zeros(matrix(1), matrix(2), 3,3);
M_error= zeros(matrix(1:2));
sqErr= zeros(matrix(1:2));

%% tensor calculation
Y1 =singleSliceDWI_set;
for x_coor =1: matrix(2)
    for y_coor =1: matrix(1)
        Y  = [];	% clear
        if mask(y_coor, x_coor) ~= 0;
            tempY = Y1(y_coor,x_coor,:);
            tempY(isnan(tempY)) = 0;
            tempY(tempY <= 0) = 20.0222; % warum 20.0222?
            y=squeeze(log(tempY )); 
            Y = [Y y'];
try
            aniso_a= iAnisoX*y;
catch
    ea_error('Wrong data / dti.bval/.bvec files specified.');
end
            YY = anisoX * aniso_a;
            M_error(y_coor,x_coor) = max(abs(YY - Y')); % MaxErr
            d_tens = [aniso_a(2),aniso_a(5),aniso_a(6);...
                aniso_a(5),aniso_a(3),aniso_a(7);...
                aniso_a(6),aniso_a(7),aniso_a(4)]; % check off_diaf elem order!!!!
            [vv, dd] = eig(-d_tens);
            eigVect(y_coor,x_coor,:,:) = vv(:,:);
            difComp(y_coor,x_coor,:,:) = dd(:,:);
            sqErr(y_coor,x_coor) = mean((YY - Y').^2);
        end
    end
end

EigenVal= sum(difComp, 3);

%end of function  do_calculation


function [gradlist, b0idx] = grads(indir,bvecf,bvalf)

gradlist=load([indir,bvecf]);
if ~(size(gradlist,2)==3)
gradlist=gradlist';
end

bvals=load([indir,bvalf]);
if size(bvals,2)==1 
   bvals=bvals'; 
end

b0idx = find(bvals<=10);




gradlist = [gradlist(:,2) -gradlist(:,1) gradlist(:,3)];
gradlist = gradlist';
gradlist(:,b0idx) = 0;




