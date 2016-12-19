function [res, errStr, oArg1]= dtdstruct_calculate(varargin)
%function [res, errStr, oArg1]= dtdstruct_calculate(commandStr [, arg1[, arg2[... [, argN]]]])
%
%
%   commandStr: {'generic'}
%   arg[1..end] the parameter specified by the command string
%
% return values:
%   res: the modified dtdStruct or if an error occured the old one
%   errStr: a message string, specifying the error which occured
% 
% Description of the different commands: 
%
%   'generic'   arg1: diffusion weighted raw data as 4D mrStruct. In the 4th
%                     dimension the diffusion encoding gradient have to vary
%               arg2: direction of the diffusion encoding gradients as
%               bTenso [Nx3x3], directions as [Nx3] or as mat or text file
%               ...
%               arg3: b-value
% 
%
% Bjoern W. Kreher
% 07/06
% Kamil Il'yasov
%
% UNIX


res= [];    errStr= '';



if (nargin < 1) || ~ischar(varargin{1})
    errStr= strcat(mfilename, '(error): First parameter have to bis a string');
    return;
end

commandStr= varargin{1};

argMax= 10;
paramCell= cell(argMax, 1);

for i= 2:nargin
    paramCell{i-1}= varargin{i};
end

if strcmp(commandStr, 'generic')
     if isempty(paramCell{5}) || ~ishandle(paramCell{5}(1))
         vHd= [];
     else
        vHd= paramCell{5}(1);
     end
     raw= []; fName= '';
     
     if ~mrstruct_istype(paramCell{1}) % rawdata
         [raw, fName, errStr]= local_import_rawdata(paramCell{1}, vHd);
         if isempty(raw)
             return;
         end
     else
         raw= paramCell{1}; fName= 'dummies.mat';
     end
 
    bVal= paramCell{3}; 

    if ~isnumeric(paramCell{2}) || (size(paramCell{2}, 2) ~= 3) || (size(paramCell{2}, 3) ~= 3)% DE direction
        [bMatrix, errStr]= local_import_bMatrix(paramCell{2}, bVal, vHd);
        if isempty(bMatrix)
            return;
        end
    else
        bMatrix= paramCell{2};
    end

    if  isnumeric(paramCell{4}) && (numel(paramCell{4}) == 1)
        threshold= paramCell{4};
    else
        errStr= strcat(mfilename, '(error): There was not threshold defined');
        return
    end

    if ~isempty(paramCell{4})
        if ischar(paramCell{4})
            fName= paramCell{4}
        else
            fName= [];
        end
    end
    [res, errStr, oArg1]= local_genericDTI_calculation(raw, bMatrix, fName, threshold, vHd);
else
    errStr= strcat(mfilename, '(error): The command ''', commandStr, ''' is not imlpemented yet');
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [raw, fName, errStr]= local_import_rawdata(fName, vHd);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raw, fName, errStr]= local_import_rawdata(fName, vHd)
raw= [];  errStr= [];

if isempty(fName)
    try
        butStr=questdlg('Do you have a mrStruct already or do you have to import data first?', ...
            'Open raw data', ...
            'mrStruct','Import','Cancel','Import');
    catch
        errStr= strcat(mfilename, ':local_import_rawdata(error): Cannot open window!');
        return
    end
    
    if strcmp(butStr, 'Import')
        msgStr= sprintf('Please select the needed dicom files and create a 4D mrStruct');
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
        
        [raw, fNameCell, dirStr]= dicom_read_tool;
        
        if isempty(raw)
            errStr= strcat(mfilename, ':local_import_rawdata(error): Aborted by user!');
            return
        end
        
        if ~mrstruct_istype(raw)
            errStr= strcat(mfilename, ':local_import_rawdata(error): result structur have to be a mrStruct');
            return
        end
        
        [p, name, ext]= fileparts(fNameCell(1, :));
        fName= fullfile(dirStr, strcat(name, '_RAW.mat'));
        
        msgStr= sprintf('Save raw data as mrStruct under %s ...', fName);
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
        
        mrstruct_write(raw, fName);
        
        msgStr= sprintf('Raw data saved under %s ...', fName);
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
    elseif strcmp(butStr, 'mrStruct')
        
        msgStr= sprintf('Open a mrStruct containing the raw data', fName);
        if isempty(vHd)
            disp(msgStr);
        else
            set(vHd, 'String', msgStr);
            drawnow;
        end
        
        [raw, fName]= mrstruct_read;
        
        if isempty(raw)
            errStr= strcat(mfilename, ':local_import_rawdata(error): No valid mrStruct');
            return
        end
    else
        errStr= strcat(mfilename, ':local_import_rawdata: Aborted by user');
        return   
    end
elseif mrstruct_istype(fName)
    raw= fName;
    if isempty(raw.patient) || ~ischar(raw.patient)
        fName= fullfile(pwd, 'dummy_RAW.mat');
    else
        fName= raw.patient;
        idx= find(fName, ' ');
        fName(idx)= '_';
        allLetter= 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_ ';  
        [X, Y]= meshgrid(fName, allLetter);                                              
        idx= sum(double(X==Y), 1) == 0;                                                  
        fName(idx)= [];
        fName= fullfile(pwd, stcat(fName, '_RAW.mat'));
    end
    
    msgStr= sprintf('Save raw data as mrStruct under %s ...', fName);
    if isempyt(vHd)
        disp(msgStr);
    else
        set(vHd, 'String', msgStr);
        drawnow;
    end
    
    mrstruct_write(raw, fName);
    
    msgStr= sprintf('Raw data saved under %s ...', fName);
    if isempyt(vHd)
        disp(msgStr);
    else
        set(vHd, 'String', msgStr);
        drawnow;
    end
end

%
%
%  START:
%       [bMatrix, errStr]= local_import_bMatrix(fName, bVal, vHd);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bMatrix, errStr]= local_import_bMatrix(fName, bVal, vHd)
bMatrix= []; errStr= [];


if ischar(fName)
    if ~exist(fName, 'file')
        errStr= strcat(mfilename, ':local_import_bMatrix: File does not exist');
        return   
    end
    
    try
        isText= 1;
        fid = fopen(fName);
        dXYZ = textscan(fid, '%n%n%n');
        fclose(fid);
        [dX, dY, dZ]= deal(dXYZ{:});        
    catch
        isText= 0;
    end
    
    if isText == 0
        try
            tmp= open(fName);
        catch
            errStr= strcat(mfilename, ':local_import_bMatrix: Could not open b-tensor file');
            return   
        end
        if isempty(tmp)
            errStr= strcat(mfilename, ':local_import_bMatrix: Could not open b-tensor file');
            return   
        end
        
        tmp= struct2cell(tmp);
        data= squeeze(tmp{1});
    else
        data= [dX dY dZ];  % bei text files nur richtungen
    end
else
    data= fName;
end

dSize= size(data);

if (length(dSize) == 3) && (dSize(1) == 3) && (dSize(2) == 3) % b-tensor schon da!!!
    bMatrix= data;
    return;
elseif (length(dSize) == 2) && (dSize(2) == 3)               % berechne b-tensor
%    data(:,2) = -data(:,2); % ika 030114  % 31.10.03 - why I have done that????

    lenAy= sqrt(sum(data.^2, 2));
    tmpAy= lenAy; 
    tmpAy(lenAy == 0)= 1;
    dirVc= data./(tmpAy*ones(1, 3));
    bMatrix =zeros(3, 3, size(dirVc, 1));
    
    if  ~isnumeric(bVal) || (numel(bVal) ~= 1)
        errStr= strcat(mfilename, '(error): Given only directions, you have to define a b-value');
        return
    end
    
    for i = 1:size(dirVc, 1)
        tmp=dirVc(i,:)'*dirVc(i,:);
        if (lenAy(i) > 0) && (sum(sum(abs(tmp) )) > 0) % ika 041118
            bMatrix(:,:,i)= (bVal/trace(tmp))*tmp;
        end
    end
else
    errStr= strcat(mfilename, ':local_import_bMatrix: Invalide data format');
    return   
end

    

%
%
%  START:
%       [res, errStr]= local_genericDTI_calculation(raw, bTensor, fName, vHd);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, dtd]= local_genericDTI_calculation(raw, bTensor, fName, threshold, vHd)
res= []; errStr= [];

if ~mrstruct_istype(raw)
    errStr= strcat(mfilename, ':local_genericDTI_calculation(error): raw have to be a valid mrStruct');
    return   
end

dwNo= size(bTensor, 3);
if (size(raw.dataAy, 4) ~= dwNo)
    errStr= strcat(mfilename, ':local_genericDTI_calculation(error): raw data are not fitting to bTensor');
    return   
end

if ~isempty(fName)
    [p, name, ext]= fileparts(fName);
    idx= find(name == '_');
    if isempty(idx)
        fName= fullfile(p, strcat(name, '_DTD.mat'));
    else
        fName= fullfile(p, strcat(name(1:idx(end)), '_DTD.mat'));
    end
end


msgStr= sprintf('Please wait ... allocate memory');
if isempty(vHd)
    %disp(msgStr);
else
    set(vHd, 'String', msgStr);
    drawnow;
end

sizeAy= mrstruct_query(raw, 'sizeAy');
b0_image_struc= mrstruct_init('volume', zeros([sizeAy(1:3)]), raw);
EigenVect= mrstruct_init('series3DEchos', zeros([sizeAy(1:3) 3 3]), raw);
EigenVal= mrstruct_init('series3D', zeros([sizeAy(1:3) 3]), raw);
error_struc= mrstruct_init('volume', zeros([sizeAy(1:3)]), raw);
[dtd, errStr]= dtdstruct_init('DTD', EigenVect, EigenVal, error_struc,'error_struc', b0_image_struc, 'b0_image_struc');

for slI= 1:sizeAy(3)
    msgStr= sprintf('Please wait ... proccessing slice %d/%d', slI, sizeAy(3));
    if isempty(vHd)
        %disp(msgStr);
    else
        set(vHd, 'String', msgStr);
        drawnow;
    end
    
    [b0Mean, eigenVal, eigVect, M_error, sqErr] = private_do_dti_calculation(bTensor, squeeze(raw.dataAy(:, :, slI, :)), threshold);
    
    dtd.b0_image_struc.dataAy(:, :, slI)= b0Mean;
    dtd.eigenVec_struc.dataAy(:, :, slI, :, :)= eigVect;
    dtd.eigenVal_struc.dataAy(:, :, slI, :)= eigenVal;
    dtd.error_struc.dataAy(:, :, slI)= M_error;
    
end
dtd.version= 'V1.1';

msgStr= sprintf('Please wait ... write file dtdStruct as %s', fName);
if isempty(vHd)
    %disp(msgStr);
else
    set(vHd, 'String', msgStr);
    drawnow;
end

if isempty(fName)
    res= 1;
else
    dtdstruct_write(dtd, fName);
    res= fName;
end

msgStr= sprintf('DTI calulation finished');
if isempty(vHd)
    %disp(msgStr);
else
    set(vHd, 'String', msgStr);
    drawnow;
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [b0Mean, EigenVal, eigVect, M_error, sqErr] = private_do_dti_calculation(B_tensor, rawSliceAy, threshold)

b0Mean= []; EigenVal= []; eigVect= []; M_error= []; sqErr= [];

nsteps =length(B_tensor);
X1 = [  reshape(B_tensor(1,1,:), 1, nsteps , 1)]; %1st diag element
X2 = [  reshape(B_tensor(2,2,:), 1, nsteps , 1)]; %2nd diag element
X3 = [  reshape(B_tensor(3,3,:), 1, nsteps , 1)]; %3rd diag element
X4 = [  reshape(B_tensor(1,2,:), 1, nsteps, 1)];
X5 = [  reshape(B_tensor(1,3,:), 1, nsteps, 1)];
X6 = [  reshape(B_tensor(2,3,:), 1, nsteps, 1)];
anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
iAnisoX=pinv(anisoX);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%   do dti calculation for the paticular slice 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bVal, dirVc]= ten2dir(B_tensor);

b0Mean= mean(rawSliceAy(:, :, bVal < 100), 3);     %%% nicht schn aber sollte erst ma funktionieren!!!!
mask= double(b0Mean > threshold);

% reserviere speicher
sizeIm= size(rawSliceAy);
%D_trace= zeros(sizeIm(1:2));
eigVect=zeros(sizeIm(1), sizeIm(2), 3,3);
difComp=zeros(sizeIm(1), sizeIm(2), 3,3);
M_error= zeros(sizeIm(1:2));
sqErr= zeros(sizeIm(1:2));


for x_coor= 1:size(mask, 2)
    for y_coor= 1:size(mask, 1)
        if (mask(y_coor, x_coor) ~= 0)
            tempY =rawSliceAy(y_coor, x_coor, :);
            itemp =find(tempY <= 0);
            if(length(itemp) >=1),
                %     'warning! signal is zero or negative',tempY(itemp),  tempY(itemp)= 20.0222; x_coor, y_coor,%  itemp,
                tempY(itemp)= 20.0222; 
            end
%            y=squeeze( sum( sum( log(tempY ), 2 ),1 ) );  % 2:5 - from 2 to 5th image
            y=squeeze(log(tempY ));  % 2:5 - from 2 to 5th image
%            anisoX=[ones(size(X1')), X1', X2', X3',2*X4', 2*X5', 2*X6'];
%            aniso_a=anisoX\y;
            aniso_a= iAnisoX*y;
            YY = anisoX*aniso_a;
%            MaxErr = max(abs(YY- y));
            
            d_tens= [aniso_a(2), aniso_a(5) ,aniso_a(6) ;aniso_a(5) ,aniso_a(3) ,aniso_a(7); aniso_a(6) ,aniso_a(7) ,aniso_a(4) ] ; % check off_diaf elem order!!!!
            [vv, dd] = eig(-d_tens) ;
            %'vv -eigen-vector, dd - eigen-value';
            eigVect(y_coor, x_coor , :,:) =vv(:,:);
            difComp(y_coor, x_coor , :,:) =dd(:,:);
            % D_raw(y_coor, x_coor,:,: ) = -d_tens;  % ika9041130
%            M_error(y_coor, x_coor ) = MaxErr;
%            sqErr(y_coor, x_coor ) = mean((YY- y).^2);
        end %if
    end % for
end % for

EigenVal= sum(difComp, 3); 
