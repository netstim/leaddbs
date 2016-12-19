% calculate_dti: calculate DTDstruct from DTI data file (Susanne Schnell)
%
%       [res,errStr] = calculate_DTI(FileName,thresholdAbs,SaveSlice,DEscheme,bfactor,vHd);
%       function saves an DTDstruct Matlab file 
%
%       example: [res,errStr] = calculate_dti('E:\DTI\Test_Diffusion\TestDebug\TestDebug_mrstruct.mat',40,0,DEscheme,1000,[]);
%
%       - needs the 'FileName' without '_raw.bin'or '_info.mat' ending saved or mrstruct-name with '.mat' ending!
%       - function 'read_DTI_Dicom' or 'read_DTI_Bruker' or 'dicom_read_tool' should have been run before
%       -> mrstruct has to be like this: [matrix(1),matrix(2),DiffEncoding,slices]
%       - SaveSlice: if equals 1 each slice is saved as mrstruct seperately, else no slice saved
%       - DEscheme needs only to be filled if '_info.mat' does not exist or if there are no Diffusion Directions saved in '_info.mat',
%       otherwise it should be empty ([])
%       - bfactor (similar to DEscheme)
%       - needs handle vHd for message display in dti_calculation_gui
%            
%       function uses tho following functions:
%           dw_data_admin
%           mrstruct_read
%           mrstruct_init
%           dtdstruct_write
%           ...
%       
%       Please take care:
%       - not tested yet
%
% code and knowledge partly copied from Kamil Ilyasov (ReadDicomFilesDTI_MZv5.m, ...)
%__________________________________________________________________________
%		Susanne Schnell
%       Department of Diagnostic Radiology, Medical Physics
%       University Hospital Freiburg
%       Hugstetter Strasse 55
%       79106 Freiburg, Germany
%		08/2007
%
%-------------------- development -------------------------------------
%   AUTHORS:
%           Susanne Schnell (SuS)
%           Kamil Il'yasov (ika)
%           Marco Reisert (mrc)
%%

function [res, msg] = calculate_dti(FileName,thresholdAbs,SaveSlice,DE_scheme,bfactor,nob0s,vHd)


%% check for input variables
if nargin > 7 || nargin < 1 
    msg= 'Error in dicom_read_singlefile: wrong number of necessary input arguments.';
    res = 0;
    return;   
end

if nargin == 1
    vHd = [];
    SaveSlice = 0;
    DE_scheme = [];
    bfactor = [];
    nob0s = [];
    thresholdAbs = 40;
elseif nargin == 2
    vHd = [];
    SaveSlice = 0;
    DE_scheme = [];
    bfactor = [];
    nob0s = [];
elseif nargin == 3
    DE_scheme = [];
    bfactor = [];
    nob0s = [];
    vHd = [];
elseif nargin == 4
    bfactor = [];
    nob0s = [];
    vHd = [];
elseif nargin == 5
    nob0s = [];
    vHd = [];
elseif nargin == 6
    vHd = [];
end

if ~ishandle(vHd)
    vHd= [];
end

%% initiate variables
EigenVal1= [];
eigVect1= [];
M_error= [];
sqErr= [];

msg= sprintf('Retrieving information ...');
if isempty(vHd)
    disp(msg);
else
    set(vHd, 'String', msg);
    drawnow;
end


%% read necessary information from file
[pathstr, name, ext] = fileparts(FileName);
if strcmp(ext,'.mat')
    dti = mrstruct_read(FileName); % mrstruct has to be like this: [matrix(1),matrix(2),slices,DiffEncoding]
else
    [res, err] = dw_data_admin('open',FileName);
    if res == 0
        error(err)
        return;
    end
end


%%   retrieve DEdir information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
% find position of b=0 scans
if strcmp(ext,'.mat')
    % check if information is missing!
    if  isempty(bfactor)
        prompt={'Enter the b-value:'};
        name='b-value and number of b0s';
        numlines=1;
        defaultanswer={'1000'};
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        answer=inputdlg(prompt,name,numlines,defaultanswer,options);
        bfactor = str2num(answer{1});
    end
% 
%     if isempty(bfactor) && isempty(nob0s)
%         prompt={'Enter the b-value:','Enter amount ob b = 0 scans:'};
%         name='b-value and number of b0s';
%         numlines=1;
%         defaultanswer={'1000','1'};
%         options.Resize='on';
%         options.WindowStyle='normal';
%         options.Interpreter='tex';
%         answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%     elseif isempty(nob0s) && ~isempty(bfactor)
%         prompt={'Enter the b-value:'};
%         name='Enter amount ob b = 0 scans:';
%         numlines=1;
%         defaultanswer={'1'};
%         options.Resize='on';
%         options.WindowStyle='normal';
%         options.Interpreter='tex';
%         answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%     elseif ~isempty(nob0s) && isempty(bfactor)
%         prompt={'Enter the b-value:'};
%         name='b-value and number of b0s';
%         numlines=1;
%         defaultanswer={'1000'};
%         options.Resize='on';
%         options.WindowStyle='normal';
%         options.Interpreter='tex';
%         answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%     end
    if isempty(DE_scheme)
        [filename, pathname,FilterIndex] = uigetfile({'*.m';'*.txt';'*.mat'},'Select file with Diffusion Encoding Scheme uncluding b0 images, has to be in image space!');
        diffDirFile = fullfile(pathname,filename);
        if FilterIndex == 0
            if isempty(vHd)
                disp('User pressed cancel, DTI calculation cannot be performed.');
            else
                set(vHd, 'String', 'Error in calculate_DTI.m: User pressed cancel, DTI calculation cannot be performed.');
                res = 0;
                drawnow;
            end
            return;
        end
        while ~exist(diffDirFile,'file')
            [filename, pathname,FilterIndex] = uigetfile({'*.m';'*.txt';'*.mat'},'File does not exist, chose again!');
            diffDirFile = fullfile(pathname,filename);
        end
        if FilterIndex == 1
            magic_Str =strcat('DE_scheme=',filename(1:end-2),';');
            eval(magic_Str)
        elseif FilterIndex == 2 | FilterIndex == 3
            DE_scheme = load(diffDirFile);        
        end
        % check if b-matrix or just listed directions! 
        if size(DE_scheme,1) == 3 && size(DE_scheme,2) == 3
            ind = find(sum(sum(DE_scheme)) == 0);
        else
            ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
        end
    else
        if isstruct(DE_scheme)
           DE_scheme = struct2cell(DE_scheme);
           DE_scheme = DE_scheme{1};
        end
        % find b0 weightings
        if size(DE_scheme,1) == 3 && size(DE_scheme,2) == 3
            ind = find(sum(sum(DE_scheme)) == 0);
        else
            ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
            % if not in DE scheme then put in b0 weightings
            if isempty(ind) && ~isempty(nob0s)
                DE_scheme = [zeros(nob0s,3);DE_scheme];
                ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
            elseif isempty(ind) && isempty(nob0s)
                error('The number of b = 0 scans is not provided, the tensor calculation cannot be started.')
                return;
            end
        end
    end
else
    if isempty(DE_scheme)
        DE_scheme = dw_data_admin('get_diffDirs');
        if ~isempty(DE_scheme)
            DE_scheme = squeeze(DE_scheme);
            ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
        else
            DE_scheme = dw_data_admin('get_bmatrix');
            ind = find(sum(sum(abs(DE_scheme))) < 300);
        end
    else
        if size(DE_scheme,1) == 3 && size(DE_scheme,2) == 3
            ind = find(sum(sum(DE_scheme)) == 0);
            if isempty(ind)
                DE_scheme = [zeros(nob0s,3,3);DE_scheme];
                ind = find(sum(sum(DE_scheme)) == 0);
            end
        else
            ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
            if isempty(ind)
                DE_scheme = [zeros(nob0s,3);DE_scheme];
                ind = find(DE_scheme(:,1) == 0 & DE_scheme(:,2) == 0 & DE_scheme(:,3) == 0);
            end
        end
    end
end

%% retrieve other information
if strcmp(ext,'.mat') %if data in mrstruct
    user = dti.user;
    TE = dti.te;
    TR = dti.tr;
    TI = dti.ti;
    edges = dti.edges;
    patient = dti.patient;
    vox = dti.vox;
    orient = dti.orient;
    matrix = size(dti.dataAy);
    b0_image = zeros(matrix(1),matrix(2), matrix(3));
    slices = matrix(3);
else % else data in binary, then info in _info.mat
    % get information from info file
    user = dw_data_admin('get_user');
    TE = dw_data_admin('get_TE');
    TR = dw_data_admin('get_TR');
    TI = dw_data_admin('get_TI');
    edges = dw_data_admin('get_edges');
    patient = dw_data_admin('get_patient');
    orient = dw_data_admin('get_orient');
    vox = dw_data_admin('get_vox');
    slices = user.SliceNo;
    b0_image = zeros(user.matrix(1),user.matrix(2),slices);
    if isempty(bfactor)
        bfactor = dw_data_admin('get_bVal');
        if isempty(bfactor)
            prompt={'Enter the b-value:'};
            name='b-value and number of b0s';
            numlines=1;
            defaultanswer={'1000'};
            options.Resize='on';
            options.WindowStyle='normal';
            options.Interpreter='tex';
            answer=inputdlg(prompt,name,numlines,defaultanswer,options);
        end
    end
end

%% calculate DTs for each slice    
for Slice = 1 : slices
    % read from file
    if strcmp(ext,'.mat')
        singleSliceDWI_set = dti.dataAy(:,:,Slice,:);
        %singleSliceDWI_set = dti.dataAy(:,:,:,Slice);
    else
        if size(DE_scheme,1) == 3 && size(DE_scheme,2) == 3
            [singleSliceDWI_set, msg] = dw_data_admin('getData', Slice, 1:size(DE_scheme,3));
        else
            [singleSliceDWI_set, msg] = dw_data_admin('getData', Slice, 1:size(DE_scheme,1));
        end
        if isempty(singleSliceDWI_set) && ~isempty(msg)
            disp(msg);
            return;
        end
    end
    singleSliceDWI_set = squeeze(singleSliceDWI_set);
    b0_image(:,:,Slice) = mean(singleSliceDWI_set(:,:,ind),3);
    if not(isfield(user,'acquNo')),
        if length(size(DE_scheme)) == 3,
            user.acquNo = size(DE_scheme,3);
        else            
            user.acquNo = size(DE_scheme,1);
        end;
    end;
    if user.acquNo > 3
        % calculate DT and DKI
        if strcmp(ext,'.mat')
            [EigenVal, eigVect, M_error, sqErr, kurtP] = DWI_calculation(DE_scheme,singleSliceDWI_set,matrix(1:2),thresholdAbs);
        else
            [EigenVal, eigVect, M_error, sqErr, kurtP] = DWI_calculation(DE_scheme,singleSliceDWI_set,user.matrix,thresholdAbs);
        end
        
        
        if isempty(EigenVal1)
            siz= size(EigenVal);
            EigenVal1= zeros(siz(1), siz(2), slices, 3);
            eigVect1= zeros(siz(1), siz(2), slices, 3, 3);
            Max_error= zeros(siz(1), siz(2), slices);
            sq_Error= zeros(siz(1), siz(2), slices);
            if not(isempty(kurtP)),
                kurtParams = zeros(siz(1), siz(2), slices,7);
            end;
        end
    end
 
    singleSliceDWI_set(:,:,ind) = [];
    if bitget(SaveSlice,1) == 1
        filename = strcat(FileName,'_slice',num2str(Slice),'.mat');
        DTI = mrstruct_init('series2D',singleSliceDWI_set);
        DTI.tr = TR;
        DTI.tr = TE;
        DTI.ti = TI;
        DTI.patient = patient;
        DTI.edges = edges;
        DTI.user = user;
        mrstruct_write(DTI,filename);
    end
    if user.acquNo > 3
        if Slice == 1
            mean_DTI = mean(singleSliceDWI_set,3);
        else
            mean_DTI(:,:,end+1) = mean(singleSliceDWI_set,3);
        end
        EigenVal1(:,:,Slice,:)= EigenVal; %ika2k1206
        eigVect1(:,: ,Slice, :,:) = eigVect;
        Max_error(:,:,Slice) = M_error;
        sq_Error(:,:,Slice)= sqErr; % square error is not used in DTDstruct, is this necessary?
        if not(isempty(kurtP)),
            kurtParams(:,:,Slice,:) = kurtP;
        end;
    end

    res = toc;

    if user.acquNo > 3
        msgStr= sprintf('Finished tensor calculation of slice %d/%d after %d min %s sec.', Slice, slices,floor(res/60), num2str(round((res - 60*floor(res/60))*10)/10));
    else
        dataAy = singleSliceDWI_set./repmat(squeeze(b0_image(:,:,Slice)),[1 1 size( singleSliceDWI_set,3)]);
        dataAy(isnan(dataAy(:))) = 0;
        meanADCnorm(:,:,Slice) = mean(dataAy,3);
        dataAy = singleSliceDWI_set;
        meanADC(:,:,Slice) = mean(dataAy,3);
        msgStr= sprintf('Finished saving meanADC and trace of slice %d/%d after %d min %s sec.', Slice, slices,floor(res/60), num2str(round((res - 60*floor(res/60))*10)/10));
    end

    if isempty(vHd)
        disp(msgStr);
    else
        set(vHd, 'String', msgStr);
        drawnow;
    end
end




%% build  raw-slices mrstruct 
if bitget(SaveSlice,2) == 1,
   
    if strcmp(ext,'.mat')
        DWI_set = dti.dataAy(:,:,:,:);     
    else
        if size(DE_scheme,1) == 3 && size(DE_scheme,2) == 3
            [DWI_set, msg] = dw_data_admin('getData',1:slices, 1:size(DE_scheme,3));
        else
            [DWI_set, msg] = dw_data_admin('getData',1:slices, 1:size(DE_scheme,1));
        end
        if isempty(DWI_set) && ~isempty(msg)
            disp(msg);
            return;
        end
    end
    
    mrstruct_HARDI = mrstruct_init('series3D',DWI_set);

    mrstruct_HARDI.vox = vox;
    mrstruct_HARDI.edges = edges;
    mrstruct_HARDI.orient = orient;
    mrstruct_HARDI.patient = patient;
    mrstruct_HARDI.te = TE;
    mrstruct_HARDI.tr = TR;
    mrstruct_HARDI.ti = TI;

    if length(size(DE_scheme)) == 3
        mrstruct_HARDI.user.bTensor = DE_scheme;
    else
        mrstruct_HARDI.user.bDir = DE_scheme';
        for k = 1:size(DE_scheme,1),
            mrstruct_HARDI.user.bTensor(:,:,k) = DE_scheme(k,:)'*DE_scheme(k,:); %*bfactor; bfactor already in DE_scheme included
        end;
    end
    mrstruct_HARDI.user.bfactor = squeeze(mrstruct_HARDI.user.bTensor(1,1,:)+mrstruct_HARDI.user.bTensor(2,2,:)+mrstruct_HARDI.user.bTensor(3,3,:));
    

end;




if strcmp(ext,'.mat')
    clear dti;
else
    dw_data_admin('close');
end



%% save as dtd_struct
b0_image_struc= mrstruct_init('volume', b0_image);
b0_image_struc.tr = TR;
b0_image_struc.te = TE;
b0_image_struc.ti = TI;
b0_image_struc.patient = patient;
b0_image_struc.orient = orient;
b0_image_struc.edges = edges;
b0_image_struc.user = user;
b0_image_struc.user.dti_calc_Ver = 'SuS_Version_1.0';
b0_image_struc.vox = vox;
if user.acquNo > 3
    error_image_struc = mrstruct_init('volume',sq_Error,b0_image_struc);
    meanDWI_image_struc = mrstruct_init('volume',mean_DTI,b0_image_struc);
    EigenVect = mrstruct_init('series3DEchos',eigVect1,b0_image_struc);
    EigenVal = mrstruct_init('series3D',EigenVal1,b0_image_struc);
    
    if not(isempty(kurtP)),
    
        K_axial = mrstruct_init('volume',kurtParams(:,:,:,1),b0_image_struc);
        K_orth = mrstruct_init('volume',kurtParams(:,:,:,2),b0_image_struc);
        D_orth = mrstruct_init('volume',kurtParams(:,:,:,3),b0_image_struc);
        D_para_int = mrstruct_init('volume',kurtParams(:,:,:,4),b0_image_struc);
        D_para_ext = mrstruct_init('volume',kurtParams(:,:,:,5),b0_image_struc);
        D_orth_ext = mrstruct_init('volume',kurtParams(:,:,:,6),b0_image_struc);
        volfrac = mrstruct_init('volume',kurtParams(:,:,:,7),b0_image_struc);

        [dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, b0_image_struc, 'b0_image_struc',meanDWI_image_struc,'meanDWI_image_struc',...
                                            K_axial,'K_axial', K_orth ,'K_orth_maxK', D_orth, 'D_orth',D_para_int ,'D_para_int' ,D_para_ext,'D_para_ext', ...
                                             D_orth_ext , 'D_orth_ext' , volfrac,'volfrac_int');
    else   
        [dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, b0_image_struc, 'b0_image_struc',meanDWI_image_struc,'meanDWI_image_struc',error_image_struc,'error_struc');
    end;
end

if strcmp(ext,'.mat')
    FileName = strcat(FileName(1:end-4),'_DTD.mat');
else
    FileName = strcat(FileName,'_DTD.mat');
end

if user.acquNo > 3
    % sort eigenvals and eigenvects
    dtd = dtdstruct_modify(dtd,'sortEigvec');

    % write dtd
    [res, errStr] = dtdstruct_write(dtd, FileName);
    if isempty(errStr)
        msgStr= sprintf('DTI calculations done, write file dtdStruct as %s', FileName);
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
    else
        msgStr= sprintf('error in calculate_dti: DTI calculations done, BUT could not write DTD');
        res = 0;
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
    end
    
    % write HARDI
    if bitget(SaveSlice,2) == 1,
        mrstruct_write(mrstruct_HARDI,[FileName(1:end-7),'HARDI.mat']);
    end;
else
    % write mrstructs
    meanADCnorm = mrstruct_init('volume',meanADCnorm,b0_image_struc);
    meanADC = mrstruct_init('volume',meanADC,b0_image_struc);
    mrstruct_write(b0_image_struc,[FileName(1:end-7),'b0image.mat']);
    mrstruct_write(meanADC,[FileName(1:end-7),'meanADC.mat']);
    mrstruct_write(meanADCnorm,[FileName(1:end-7),'meanADCnorm.mat']);
end



%% finish...
msg= sprintf('DTI calculation finished.');
if isempty(vHd)
    disp(msg);
else
    set(vHd, 'String', msg);
    drawnow;
end

res = 1;
%***************************end of main ***********************************


%% local function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        tmp = DEscheme(i,:)' * DEscheme(i,:);  % *bvalue;
%        if sum(sum(abs(tmp) )) > 0 % ika 041118    bvalue is now included in vectorial DEscheme vias sqrt(b)        
        B_tensor(:,:,i)=  tmp;
%        end
    end
end

U = getDiffDirTrans;

for k = 1:size(B_tensor,3),
    B_tensor(:,:,k) =  U*B_tensor(:,:,k)*U';
end;



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
            aniso_a= iAnisoX*y;
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


