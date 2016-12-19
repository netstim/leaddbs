% creation of a binary data file with DWI data and an info .mat file (Bj�rn Kreher and Susanne Schnell)
%
%function [res, errStr, res2]= dw_data_admin(commandStr, [arg1[, arg2[, ...])
%
%   Method to create, write and read diffusion weighte data for fast access
%
% create: matSize, sliceNo, acquNo, fName
% open: fName
% close: 
% free:
% putData: sliceNo, acquNo, dataAy
% addData: sliceNo, acquNo, dataAy
% getData: sliceNo, acquNo
% 
% set_diffDirs: bDirVc
% get_diffDirs: 
% set_bmatrix: bmatr
% get_bmatrix: 
% set_bVal: bval (several b values possible)
% get_bVal: 
% set_edges: edgesAy
% get_edges:
% get_TR: 
% set_TR: tr
% get_TE: 
% set_TE: te
% get_TI: 
% set_TI: ti
% get_orient:
% set_orient: orient
% get_patient:
% set_patient: patient
% get_vox: 
% set_vox: vox
% get_user: 
% set_user: user
% 
%
%
% Bjoern W. Kreher
% 11/02
%
% UNIX

function [res, errStr, oArg1, oArg2]= dw_data_admin(varargin)
res= [];    errStr= ''; oArg1= [];    oArg2= [];   
maxArg= 10;

persistent rawFileHd_per fName_per dataDimension_per;


%% parameter-check and -moving
if length(varargin) < 1
    errStr= sprintf('%s(varagrin): There have to be at least two parameters', mfilename);
    return;
end
if ~ischar(varargin{1})
    errStr= sprintf('%s(varagrin): Command have to be a string', mfilename);
    return;
else
    commandStr= varargin{1};
end

param= cell(maxArg,1);
for i= 1:maxArg
    if length(varargin) < (i + 1)
        param{i}= [];
    else
        param{i}= varargin{i + 1};
    end
end

%%%%%
%%%%% begin of command switching part
if strcmp(commandStr, 'create') %% initialisiert und kreiert datenstruktur
    if ~isempty(rawFileHd_per)
         butStr=questdlg('dw_data_admin currently processes a file!', 'dw_data_admin', ...
             'Continue','Abort','Abort');
         if strcmp(butStr, 'Abort')
             errStr= sprintf('%s(create): Aborted by user', mfilename);
             return
         end
    end
    if isnumeric(param{1}) && isscalar(param{2}) && isscalar(param{3}) && ischar(param{4})
        matSize= param{1};        sliceNo= param{2};    acquNo= param{3};
        fName= param{4}; %optional
    else
        errStr= sprintf('%s(create): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr, rawFileHd_per, fName_per, dataDimension_per]= local_createFile(fName, matSize, sliceNo, acquNo);
%%%
elseif strcmp(commandStr, 'open') %% open a file which already exists
    if ~isempty(rawFileHd_per)
         butStr=questdlg('dw_data_admin currently processes a file!', 'dw_data_admin', ...
             'Continue','Abort','Abort');
         if strcmp(butStr, 'Abort')
             errStr= sprintf('%s(open): Aborted by user', mfilename);
             return
         end
    end
               
    if ischar(param{1})
        fName= param{1};
    else
        errStr= sprintf('%s(open): Wrong parameter format', mfilename);
        return        
    end        
    [res, errStr, rawFileHd_per, fName_per, dataDimension_per]= local_openFile(fName);
%%%
elseif strcmp(commandStr, 'close') %% closes the file structure which is currently open    
    [res, errStr, rawFileHd_per, fName_per, dataDimension_per]= local_closeFile(rawFileHd_per, fName_per, dataDimension_per);
%%%
elseif strcmp(commandStr, 'free') %% closes the file structure which is currently open
    [res, errStr, rawFileHd_per, fName_per, dataDimension_per]= local_freeFile(rawFileHd_per, fName_per, dataDimension_per);

elseif strcmp(commandStr, 'putData') %% closes the file structure which is currently open
    if isnumeric(param{1}) && isnumeric(param{2}) && isnumeric(param{3})
        sliceNo= param{1};  acquNo= param{2};   dataAy= param{3};
    else
        errStr= sprintf('%s(putData): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_putData(rawFileHd_per, fName_per, dataDimension_per, sliceNo, acquNo, dataAy, 'put');
elseif strcmp(commandStr, 'addData') %% closes the file structure which is currently open
    if isnumeric(param{1}) && isnumeric(param{2}) && isnumeric(param{3})
        sliceNo= param{1};  acquNo= param{2};   dataAy= param{3};
    else
        errStr= sprintf('%s(addData): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_putData(rawFileHd_per, fName_per, dataDimension_per, sliceNo, acquNo, dataAy, 'add');
elseif strcmp(commandStr, 'getData') %% closes the file structure which is currently open
    if isnumeric(param{1}) && isnumeric(param{2})
        sliceNo= param{1};  acquNo= param{2};
    else
        errStr= sprintf('%s(getData): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_getData(rawFileHd_per, fName_per, dataDimension_per, sliceNo, acquNo);    
%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%
elseif strcmp(commandStr, 'set_diffDirs') %% writes the diffusion encoding directions into the info file
    if isnumeric(param{1})
        bDirVc= param{1};
    else
        errStr= sprintf('%s(set_diffDirs): Wrong parameter format', mfilename);
        return        
    end
     [res, errStr]= local_setDiffDir(fName_per, bDirVc);
elseif strcmp(commandStr, 'get_diffDirs') %% reads the diffusion encoding directions from the info file
    [res, errStr]= local_getDiffDir(fName_per);
    
elseif strcmp(commandStr, 'set_bmatrix') %% writes the diffusion encoding directions into the info file
    if isnumeric(param{1})
        bmatr= param{1};
    else
        errStr= sprintf('%s(set_bmatrix): Wrong parameter format', mfilename);
        return
    end
    [res, errStr]= local_setbmatrix(fName_per, bmatr);
elseif strcmp(commandStr, 'get_bmatrix') %% reads the diffusion encoding directions from the info file
    [res, errStr]= local_getbmatrix(fName_per);

elseif strcmp(commandStr, 'set_bVal') %% writes the b-Value
    if isnumeric(param{1})
        bDirVc= param{1};
    else
        errStr= sprintf('%s(set_bVal): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_setBVal(fName_per, bDirVc);
elseif strcmp(commandStr, 'get_bVal') %% reads the b-Value
    [res, errStr]= local_getBVal(fName_per);    
elseif strcmp(commandStr, 'set_edges') %% writes edges
    if isnumeric(param{1})% & isnumeric(param{2}) % SuS: weitere Dim f�r edges
        edgesAy= param{1};
       % acquNo = param{2};
    else
        errStr= sprintf('%s(set_edges): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_setEdges(fName_per, edgesAy); %, acquNo); % SuS weitere Dim f�r edges von allen Bildern
elseif strcmp(commandStr, 'get_edges') %% reads edges
    [res, errStr]= local_getEdges(fName_per);    
elseif strcmp(commandStr, 'get_TR') %% reads TR
    [res, errStr]= local_getTR(fName_per);    
elseif strcmp(commandStr, 'set_TR') %% writes TR
    if isscalar(param{1})
        tr= param{1};
    else
        errStr= sprintf('%s(set_TR): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_setTR(fName_per, tr);
elseif strcmp(commandStr, 'get_TE') %% reads TE
    [res, errStr]= local_getTE(fName_per);    
elseif strcmp(commandStr, 'set_TE') %% writes TE
    if isscalar(param{1})
        te= param{1};
    else
        errStr= sprintf('%s(set_TE): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_setTE(fName_per, te);
elseif strcmp(commandStr, 'get_TI') %% reads TI
    [res, errStr]= local_getTI(fName_per);    
elseif strcmp(commandStr, 'set_TI') %% writes TI
    if isscalar(param{1})
        ti= param{1};
    else
        errStr= sprintf('%s(set_TI): Wrong parameter format', mfilename);
        return        
    end
    [res, errStr]= local_setTI(fName_per, ti);
elseif strcmp(commandStr, 'get_orient') %% reads TI
    [res, errStr]= local_getOrient(fName_per);    
elseif strcmp(commandStr, 'set_orient') %% writes TI
    if ischar(param{1})
        orient= param{1};
    else
        errStr= sprintf('%s(set_orient): Parameter is not a string', mfilename);
        return        
    end
    [res, errStr]= local_setOrient(fName_per, orient);
elseif strcmp(commandStr, 'get_patient') %% reads TI
    [res, errStr]= local_getPatient(fName_per);    
elseif strcmp(commandStr, 'set_patient') %% writes TI
    if ischar(param{1})
        patient= param{1};
    else
        patient = num2str(param{1});
        if isempty(patient)
            patient = 'no patient name';
        end
%         errStr= sprintf('%s(set_patient): Parameter is not a string', mfilename);
%         return        
    end
    [res, errStr]= local_setPatient(fName_per, patient);
elseif strcmp(commandStr, 'get_vox') %% reads TI
    [res, errStr]= local_getVox(fName_per);    
elseif strcmp(commandStr, 'set_vox') %% writes TI
    if isnumeric(param{1})
        vox= param{1};
    else
        errStr= sprintf('%s(set_vox): Parameter is not a string', mfilename);
        return        
    end
    [res, errStr]= local_setVox(fName_per, vox);
elseif strcmp(commandStr, 'get_user') %% reads TI
    [res, errStr]= local_getUser(fName_per);    
elseif strcmp(commandStr, 'set_user') %% writes TI
    if isstruct(param{1})
        user= param{1};
    else
        errStr= sprintf('%s(set_user): Parameter is not a struct', mfilename);
        return        
    end
    [res, errStr]= local_setUser(fName_per, user);
else
    errStr= sprintf('%s(varagrin): command ''%s'' is not implemented yet', mfilename, commandStr);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setUser(fName, user)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.user= user;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getUser(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.user;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%

    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setVox(fName, vox)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

len= numel(vox);

if (len == 3)
    infoStrc.vox= [vox(1) vox(2) vox(3) 0];
elseif (len == 4)
    infoStrc.vox= [vox(1) vox(2) vox(3) vox(3)];
else
    errStr= sprintf('%s::local_setVox(error): wrong format of vox', mfilename);
    return;    
end
[res, errStr]= private_writeInfoStruct(fName, infoStrc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [res, errStr]= local_getVox(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.vox;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setPatient(fName, patient)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.patient= patient;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getPatient(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.patient;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setOrient(fName, orient)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.orient= orient;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getOrient(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.orient;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setEdges(fName, edges) %,acquNo) % weitere Dim f�r edges von allen Bildern
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

if  ~isequal(size(edges), [4, 4]) %,acquNo]) %SuS edges von allen Bildern
    errStr= sprintf('%s::local_setEdges(error): wrong format of edges matrix (4x4)', mfilename);
    return;    
end
infoStrc.edges= edges;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [res, errStr]= local_getEdges(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.edges;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setTR(fName, tr)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.TR= tr;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getTR(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.TR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setTE(fName, te)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.TE= te;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getTE(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.TE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setTI(fName, ti)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
infoStrc.TI= ti;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getTI(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.TI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%
    
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setBVal(fName, bVal)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

len= numel(bVal);

%if (len ~= 1) & (len ~= infoStrc.volNo)
%    errStr= sprintf('%s::local_setBVal(error): wrong format of bValue', mfilename);
%    return;    
%end
infoStrc.bValue= reshape(bVal, [len, 1]);
[res, errStr]= private_writeInfoStruct(fName, infoStrc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [res, errStr]= local_getBVal(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.bValue;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setDiffDir(fName, bDirVc)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

sizeAy= size(bDirVc);

if ~isequal(sizeAy, [infoStrc.volNo 3])
    errStr= sprintf('%s::local_setDiffDir(error): wrong number of dw-directions', mfilename);
    return;    
end
infoStrc.dirDirVc= bDirVc;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getDiffDir(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.dirDirVc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_setbmatrix(fName, bmatr)
res= []; errStr= '';
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

% CHECK THIS!
sizeAy= size(bmatr);

if ~isequal(sizeAy, [3 3 infoStrc.volNo])
    errStr= sprintf('%s::local_setbmatrix(error): wrong dimensions must be a 3x3xDEdirection matrix', mfilename);
    return;    
end
infoStrc.bmatrix= bmatr;
[res, errStr]= private_writeInfoStruct(fName, infoStrc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [res, errStr]= local_getbmatrix(fName)
res= []; errStr= '';

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
res= infoStrc.bmatrix;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%

% res.matSize=                [];     %Size of the matrix [int, int]
% res.sliceNo=                [];     %Number or slices [int]
% res.volNo=                  [];     %Nummer of acquired volumes [int]
% res.bValue=                 [];     %either one b-value or array of b-values [double] or [volNo x 1]
% 
% res.edges=                  [];     %homogene transformation matrix [4x4]
% res.TR=                     [];     % tR [double]
% res.TE=                     [];     % te [double]
% res.TI=                     [];     % tI [double]
% res.patient=                '';     %name of patient
% res.orient=                 [];     %patient orientation 
% res.vox=                    [];     %voxel size
% res.user=                   [];     %Structure of additional data



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [res, errStr]= local_getData(rawHd, fName, dataDim, sliceNo, acquNo)
res= [];    errStr= '';   

if isempty(rawHd)
    errStr=  sprintf('%s::local_getData(error): no file is open', mfilename);
    return;
end

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end
if ~isequal(sliceNo, round(abs(sliceNo))) || (min(sliceNo) <= 0) || (max(sliceNo) > dataDim(3))
    errStr=  sprintf('%s::local_getData(error): no valid slice number', mfilename);
    return;    
end
if ~isequal(acquNo, round(abs(acquNo))) || (min(acquNo) <= 0) || (max(acquNo) > dataDim(4))
    errStr=  sprintf('%s::local_getData(error): no valid volume number', mfilename);
    return;    
end

dataAy= zeros(dataDim(1), dataDim(2), length(sliceNo), length(acquNo));

sd= sizeofDouble;
count= 0;
for i= 1:length(sliceNo)
    for j= 1:length(acquNo)
        if infoStrc.repNo(sliceNo(i), acquNo(j)) < 1
            warning('%s::local_getData(warning): Image of slice %d and acquNo %d contains no data', mfilename,sliceNo(i), acquNo(j));
        end
        offset= dataDim(1)*dataDim(2)*(sliceNo(i) - 1) + dataDim(1)*dataDim(2)*dataDim(3)*(acquNo(j) - 1);
        
        if fseek(rawHd, offset*sd, 'bof') < 0
            errStr= sprintf('%s::local_getData(error): could not move to data pos (%d; %d)', ...
                mfilename, sliceNo(i), axquNo(j));
            return;
        end
        
        [tmpAy, no]= fread(rawHd, [dataDim(1) dataDim(2)], 'double');
        if no ~= dataDim(1)*dataDim(2)
            errStr= sprintf('%s::local_getData(error): could write data at sliceId: %d; acquId: %d', ...
                mfilename, sliceNo(i), axquNo(j));
            return;
        end
        %dataAy(:, :, i, j)= tmpAy/infoStrc.repAy(sliceNo(i), acquNo(j));
        dataAy(:, :, i, j)= tmpAy/infoStrc.repNo(sliceNo(i), acquNo(j));
        count= count + no;
    end
end

res= dataAy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= local_putData(rawHd, fName, dataDim, sliceNo, acquNo, dataAy, addFlag)
res= -1;    errStr= '';   

sizeAy= size(dataAy);
sizeAy((end+1):4)= 1;

[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

if isempty(rawHd)
    errStr=  sprintf('%s::local_putData(error): no file is open', mfilename);
    return;
end
if ~isequal(sizeAy([1 2]), dataDim([1 2]))
    errStr=  sprintf('%s::local_putData(error): wrong matrix size', mfilename);
    return;
end
if (sizeAy(3) ~= length(sliceNo)) || (sizeAy(4) ~= length(acquNo))
    errStr=  sprintf('%s::local_putData(error): wrong data dimension ', mfilename);
    return;
end

if ~isequal(sliceNo, round(abs(sliceNo))) || (min(sliceNo) <= 0) || (max(sliceNo) > dataDim(3))
    errStr=  sprintf('%s::local_putData(error): no valid slice number', mfilename);
    return;    
end
if ~isequal(acquNo, round(abs(acquNo))) || (min(acquNo) <= 0) || (max(acquNo) > dataDim(4))
    errStr=  sprintf('%s::local_putData(error): no valid volume number', mfilename);
    return;    
end

sd= sizeofDouble;
count= 0;
for i= 1:length(sliceNo)
    for j= 1:length(acquNo)
        offset= dataDim(1)*dataDim(2)*(sliceNo(i) - 1) + dataDim(1)*dataDim(2)*dataDim(3)*(acquNo(j) - 1);
        
        if fseek(rawHd, offset*sd, 'bof') < 0
            errStr= sprintf('%s::local_putData(error): could not move to data pos (%d; %d)', ...
                mfilename, sliceNo(i), acquNo(j));
            return;
        end
        
        if strcmp(addFlag, 'add') && (infoStrc.repNo(sliceNo(i), acquNo(j)) > 0) % Bei mittelungen k�nnen daten auch aufaddiert werden
            [tmpAy, no]= fread(rawHd, [dataDim(1) dataDim(2)], 'double');
            if no ~= dataDim(1)*dataDim(2)
                errStr= sprintf('%s::local_putData(error)[add]: could read data at sliceId: %d; acquId: %d', ...
                            mfilename, sliceNo(i), acquNo(j));
                return;
            end
            if fseek(rawHd, offset*sd, 'bof') < 0
                errStr= sprintf('%s::local_putData(error)[add]: could not move to data pos (%d; %d)', ...
                    mfilename, sliceNo(i), acquNo(j));
                return;
            end
            no= fwrite(rawHd, dataAy(:, :, i, j) + tmpAy, 'double');
            infoStrc.repNo(sliceNo(i), acquNo(j))= infoStrc.repNo(sliceNo(i), acquNo(j)) + 1;
        else
            no= fwrite(rawHd, dataAy(:, :, i, j), 'double');
            infoStrc.repNo(sliceNo(i), acquNo(j))= 1;
        end
        if no ~= dataDim(1)*dataDim(2)
            errStr= sprintf('%s::local_putData(error): could write data at sliceId: %d; acquId: %d', ...
                mfilename, sliceNo(i), acquNo(j));
            return;
        end
        count= count + no;
    end
end

[res, errStr]= private_writeInfoStruct(fName, infoStrc);
if res <= 0
    res= [];
    return;
end
res= count;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, rawFileHdR, fNameR, dataDimR]= local_closeFile(rawFileHd_per, fName_per, dataDimension_per)
res= -1;    errStr= '';   fNameR= ''; rawFileHdR= [];     dataDimR= [];

% check arguments
if isempty(rawFileHd_per)
    errStr=  sprintf('%s::local_closeFile(error): no file is open', mfilename);
    return;
end

try
    res= fclose(rawFileHd_per);
catch
    res= -1;
    lerr= lasterr;
    lerr(lerr == 10)= 59;
    errStr=  sprintf('%s::local_closeFile: was not able to close file (%s)', mfilename, fName, lerr);
    return 
end
res= 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, rawFileHdR, fNameR, dataDimR]= local_freeFile(rawFileHd_per, fName_per, dataDimension_per)
res= [];    errStr= '';   fNameR= ''; rawFileHdR= [];     dataDimR= [];

[res, errStr, rawFileHdR, fNameR, dataDimR]= local_closeFile(rawFileHd_per, fName_per, dataDimension_per);

if res < 0
    return;
end

try
    delete(strcat(fName_per, '_info.mat'));
    delete(strcat(fName_per, '_raw.bin'));
catch
    res= -2;
    lerr= lasterr;
    lerr(lerr == 10)= 59;
    errStr=  sprintf('%s::local_closeFile: was not able to close file (%s)', mfilename, fName, lerr);
    return 
end
res= 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, rawFileHdR, fNameR, dataDimR]= local_openFile(fName)
res= [];    errStr= '';     fNameR= ''; rawFileHdR= [];     dataDimR= [];

%%%%%%%% open infostruct
[infoStrc, errStr]= private_readInfoStruct(fName);
if isempty(infoStrc)
    return;
end

dataDim= [infoStrc.matSize(1), infoStrc.matSize(2), infoStrc.sliceNo, infoStrc.volNo];

%%%%%%% open raw_file
fName_raw= strcat(fName, '_raw.bin');
rawFileHd= fopen(fName_raw,'r+', 'ieee-le');
if rawFileHd < 0
    errStr= sprintf('%s::local_openFile(error): could not open raw file', mfilename);
    return;
end

%%%%%%% check rawfile size
fseek(rawFileHd, 0, 'eof');
count= ftell(rawFileHd);

if count ~= (prod(dataDim)*sizeofDouble)
    errStr= sprintf('%s::local_openFile(error): file has wrong length', mfilename);
    fclose(rawFileHd);
    return;
end

fNameR= fName; res= count;    rawFileHdR= rawFileHd;     dataDimR= dataDim;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, rawFileHdR, fNameR, dataDimR]= local_createFile(fName, matSize, sliceNo, acquNo)
res= [];    errStr= '';     fNameR= ''; rawFileHdR= [];     dataDimR= [];

% check arguments
if (numel(matSize) ~= 2) && isequal(matSize, abs(round(matSize)))
    errStr=  sprintf('%s::local_createFile(error): matSize has wrong format', mfilename);
    return;
end
if sliceNo ~=  abs(round(sliceNo))
    errStr=  sprintf('%s::local_createFile(error): sliceNo has wrong format', mfilename);
    return;
end
if acquNo ~= abs(round(acquNo))
    errStr=  sprintf('%s::local_createFile(error): acquNo has wrong format', mfilename);
    return;
end

dataDim= [matSize(1), matSize(2), sliceNo, acquNo];
if (prod(dataDim) <= 0)
    errStr=  sprintf('%s::local_createFile(error): data dimension is not valid', mfilename);
    return;
end

%%%%%%%% create infostruct
infoStrc= private_initInfoStruct;
infoStrc.matSize= [matSize(1), matSize(2)];
infoStrc.sliceNo= sliceNo;
infoStrc.volNo=   acquNo;
infoStrc.repNo= zeros(sliceNo, acquNo);
[res, errStr]= private_writeInfoStruct(fName, infoStrc);
if res <= 0
    res= [];
    return;
end

%%%%%%% create raw_file
fName_raw= strcat(fName, '_raw.bin');
rawFileHd= fopen(fName_raw,'w+', 'ieee-le');
if rawFileHd < 0
    errStr= sprintf('%s::local_createFile(error): could not create raw file', mfilename);
    return;
end
%%%%%%% reserver sape on disk
tmpAy= zeros(matSize(1), matSize(2), sliceNo);
count= 0;
for i= 1:acquNo
    no= fwrite(rawFileHd, tmpAy, 'double');
    if no ~= numel(tmpAy)
        errStr= sprintf('%s::local_createFile(error): could not reach end of file', mfilename);
        return;
    end
    count= count + no;
end

fNameR= fName; res= count;    rawFileHdR= rawFileHd;     dataDimR= dataDim;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function res= sizeofDouble
%%%%%% determine size of double
tmp= double(1);
wh= whos('tmp');
res= wh.bytes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_initInfoStruct
errStr= '';

res.matSize=    [];     %Size of the matrix [int, int]
res.sliceNo=    [];     %Number or slices [int]
res.volNo=      [];     %Nummer of acquired volumes [int]
res.bValue=     [];     %either one b-value or array of b-values [double] or [volNo x 1]
res.dirDirVc=   [];     %Diffusion encoding directions / b0 images are marked by vector [0 0 0] [volNo x 3]
res.edges=      [];     %homogene transformation matrix [4x4]
res.TR=         [];     % tR [double]
res.TE=         [];     % te [double]
res.TI=         [];     % tI [double]
res.patient=    '';     %name of patient
res.orient=     [];     %patient orientation 
res.vox=        [];     %voxel size
res.user=       [];     %Structure of additional data
res.repNo=      [];     %array containing the number of repetition per image (to calculate the mean) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_IsInfoStruct(infoStrc)
res= false; errStr= '';

if ~isstruct(infoStrc)
    errStr=  sprintf('%s::private_IsInfoStruct: param is no struct', mfilename);
    return;
end

membersCell= {'repNo'; 'matSize'; 'sliceNo'; 'volNo'; 'bValue'; 'dirDirVc'; 'edges'; 'TR'; 'TE'; 'TI'; 'patient'; 'orient'; 'vox'; 'user'};

for i= 1:length(membersCell)
    if ~isfield(infoStrc, membersCell{i})
        errStr=  sprintf('%s::private_IsInfoStruct: param does not contain the field ''%s''', mfilename, membersCell{i});
        return;
    end
end
res= true;
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_readInfoStruct(fName)
errStr= '';

fName= strcat(fName, '_info.mat');
try
    res= open(fName);
catch
    res= [];
    lerr= lasterr;
    lerr(lerr == 10)= 59;
    errStr=  sprintf('%s::private_readInfoStruct: was not able to open infofile ''%s'' => ''%s''', mfilename, fName, lerr);
    return 
end

[ok, errStr]= private_IsInfoStruct(res);
if ~ok
    res= [];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr]= private_writeInfoStruct(fName, infoStr)
res= false; errStr= '';

[res, errStr]= private_IsInfoStruct(infoStr);
if ~res
    return;
end

fName= strcat(fName, '_info.mat');
try
    save(fName, '-struct', 'infoStr');
catch
    res= false;
    lerr= lasterr;
    lerr(lerr == 10)= 59;
    errStr=  sprintf('%s::private_writeInfoStruct: was not able to write infofile ''%s'' => ''%s''', mfilename, fName, lerr);
    return 
end
res= true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


