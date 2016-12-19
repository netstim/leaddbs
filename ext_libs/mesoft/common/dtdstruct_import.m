function [res, errStr, oArg1, oArg2]= dtdstruct_import(varargin)
%   function [dtdStrcut, errStr, oArg1, oArg2]= dtdstruct_import(typeStr, [, op1[, op2[... [, opN]]]])
%
%   command:    {'SPM2' | '' | 
%
%   'SPM2': [transf_My, [in_fNames]]
%           transf_My: 4x4 Matrix; M*vox -> world system
%
% Bjoern W. Kreher
% 08/05
%
% UNIX


res= []; errStr= ''; oArg1= []; oArg2= [];
maxArg= 10;

verStr= strcat(mfilename, '_V0.1');

if (nargin >= 1) && ischar(varargin{1})
    typeStr= varargin{1};
else
    errStr= strcat(mfilename, ' (error): Second argument have to be of the type string');
    return
end

argCell= cell(1, maxArg);
for i= 1:(nargin - 1)
    argCell{i}= varargin{1 + i};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SPM02', typeStr) || strcmp('SPM02_SL', typeStr)
    if isempty(argCell{1}) && isempty(argCell{2}) && isempty(argCell{3}) && isempty(argCell{4})
        [argCell{1}, errStr, argCell{2}, argCell{3}, argCell{4}]= local_import_SPM2_getParam;
        if isempty(argCell{1})
            return
        end
    end
    
    if isempty(argCell{1}) || ~iscell(argCell{1}) || (length(argCell{1}) ~= 4)
        errStr= strcat(mfilename, '::SPM02 (error): First arg should contain four paths (first eigenvalue and eigenvector images)');
        return
    else
        eig1StrLst= argCell{1};        
    end
    if isempty(argCell{2}) || ~iscell(argCell{2}) || (length(argCell{2}) ~= 4)
        errStr= strcat(mfilename, '::SPM02 (error): Second arg should contain four paths (second eigenvalue and eigenvector images)');
        return
    else
        eig2StrLst= argCell{2};        
    end
    if isempty(argCell{3}) || ~iscell(argCell{3}) || (length(argCell{3}) ~= 4)
        errStr= strcat(mfilename, '::SPM02 (error): Third arg should contain four paths (third eigenvalue and eigenvector images)');
        return
    else
        eig3StrLst= argCell{3};        
    end
    if isempty(argCell{4}) && ~iscell(argCell{4})
        errStr= strcat(mfilename, '::SPM02 (error): Fourth arg should contain paths of other data modalities');
        return
    else
        otherStrLst= argCell{4};        
    end
    
    if strcmp('SPM02', typeStr)
        [res, errStr]= local_import_SPM2(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst);
        if ~isempty(res) && (nargout == 3)
            [slStrc, errStr, dirStr]= local_create_SPM2_SL(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst);        
            oArg1= {slStrc; dirStr};
        end
    else
        [res, errStr, oArg1]= local_create_SPM2_SL(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst);        
    end        
elseif strcmp('priv_SPM02_SL', typeStr)
    if isempty(argCell{1}) && ischar(argCell{1}) && (exist(argCell{1}, 'file') == 0)
        errStr= strcat(mfilename, '::SPM02_SL (error): First arg should contain the path to the softlink dtdStruct(SPM02) file');
        return
    else
        fName= argCell{1};        
    end
        
    [res, errStr]= local_import_SPM2_SL(fName);
elseif strcmp('mrStruct', typeStr)
else
    errStr= strcat(mfilename, ' (error): Data type ''', typeStr, ''' is not supported yet');
    return    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%  START:
%       [dtd, errStr]= loacal_import_SPM2(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dtd, errStr]= local_import_SPM2(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst)

dtd= []; errStr= ''; dtdStruct= []; 
dtdStruct= [];


if ischar(eig1StrLst{1}) && ischar(eig2StrLst{1}) && ischar(eig3StrLst{1})
    [eigVal_mrS, errStr]= mrstruct_import('SPM02', 'series3D', {eig1StrLst{1}; eig2StrLst{1}; eig3StrLst{1}});
    if isempty(eigVal_mrS)
        return;
    end
    idx= ~isfinite(eigVal_mrS.dataAy);
    eigVal_mrS.dataAy(idx)= 0;
else
    errStr= strcat(mfilename, 'local_import_n_SPM2 (error): Paramater have contain a cell of strings');
    return
end

 
if ischar(eig1StrLst{2}) && ischar(eig2StrLst{2}) && ischar(eig3StrLst{2}) && ...
        ischar(eig1StrLst{3}) && ischar(eig2StrLst{3}) && ischar(eig3StrLst{3}) &&...
        ischar(eig1StrLst{4}) && ischar(eig2StrLst{4}) && ischar(eig3StrLst{4})
    [eigVec_mrS, errStr]= mrstruct_import('SPM02', 'series3DEchos', ...
           {eig1StrLst{2}  eig2StrLst{2}  eig3StrLst{2}; ...
            eig1StrLst{3}  eig2StrLst{3}  eig3StrLst{3}; ...
            eig1StrLst{4}  eig2StrLst{4}  eig3StrLst{4}});
    if isempty(eigVec_mrS)
        return;
    end
    idx= ~isfinite(eigVec_mrS.dataAy);
    eigVec_mrS.dataAy(idx)= 0;
else
    errStr= strcat(mfilename, 'local_import_n_SPM2 (error): Paramater have contain a cell of strings');
    return
end

% Transform eigenvectores
sizeAy= mrstruct_query(eigVec_mrS, 'sizeAy');

afTrVc= spm_imatrix(eigVec_mrS.user.hMatrix);
transEVecMx= spm_matrix([0 0 0 afTrVc([4 5 6]) sign(afTrVc([7 8 9])) 0 0 0]);  
transEVecMx= inv(transEVecMx(1:3, 1:3));

transEVecMx= transEVecMx*diag([-1 1 1]);  %%%%%%%%%%%%%%%%%%%%%%%%%% XXX da stimmt noch was mit den Daten nicht
                                          %%%%%%%%%%%%%%%%%%%%%%%%%% XXX 05-09-02 mit Volkmar besprochen. Er wird es in SPM05 berarbeiten!
                
eVecData= reshape(permute(eigVec_mrS.dataAy, [4 5 1 2 3]), [3 3*prod(sizeAy(1:3))]);
eVecData= transEVecMx(1:3, 1:3)*eVecData;
eigVec_mrS.dataAy= permute(reshape(eVecData, [sizeAy([4 5 1 2 3])]), [3 4 5 1 2]);

% import other datas 
commandStr= '[dtd, errStr]= dtdstruct_init(eigVec_mrS, eigVal_mrS';

oMrCell= cell(length(otherStrLst), 1);
for i= 1:length(otherStrLst)
    [oMrCell{i}, errStr]= mrstruct_import('SPM02', 'volume', otherStrLst{i});
    [dummy1, fName]= fileparts(otherStrLst{i});
    commandStr= sprintf('%s, oMrCell{%d}, ''%s''',commandStr, i, str2filename(fName, 'struct'));
    idx= ~isfinite(oMrCell{i}.dataAy);
    oMrCell{i}.dataAy(idx)= 0;
    
end
commandStr= sprintf('%s);',commandStr);
eval(commandStr);

dtd.version= 'V1.1_SPM';


%
%
%  START:
%    [res, errStr]= local_create_SPM2_SL(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res, errStr, oArg1]= local_create_SPM2_SL(eig1StrLst, eig2StrLst, eig3StrLst, otherStrLst)

res= []; errStr= ''; oArg1= '';

pathSollStr= fileparts(eig1StrLst{1});
scLst= [];
fNamesCell= cat(2, reshape(eig1StrLst, [4 1]), reshape(eig2StrLst, [4 1]), reshape(eig3StrLst, [4 1]));
scLst.PRIV_FORMAT= 'SPM02';
scLst.PRIV_DATE= date;
%%%%%   length(otherStrLst);
for i= 1:3 
    %Eigenvalues
    [pathStr, fName, ext]= fileparts(fNamesCell{1, i}); 
    if ~strcmp(pathSollStr, pathStr)
        errStr= strcat(mfilename, 'local_create_SPM2_SL (error): All data files have to be in the same directory for soft link dtdStruct');
        return
    end
    scLst.(strcat('EIGVAL_', num2str(i))) = strcat(fName, ext);
    %EigenVector
    compStr= 'XYZ';
    for ii= 1:3
        [pathStr, fName, ext]= fileparts(fNamesCell{ii+1, i}); 
        if ~strcmp(pathSollStr, pathStr)
            errStr= strcat(mfilename, 'local_create_SPM2_SL (error): All data files have to be in the same directory for soft link dtdStruct');
            return
        end
        scLst.(strcat('EIGVECT_', num2str(i), '_', compStr(ii))) = strcat(fName, ext);
    end
end 
for i= 1:length(otherStrLst)
    [pathStr, fName, ext]= fileparts(otherStrLst{i}); 
    if ~strcmp(pathSollStr, pathStr)
        errStr= strcat(mfilename, 'local_create_SPM2_SL (error): All data files have to be in the same directory for soft link dtdStruct');
        return
    end
    scLst.(strcat('VOLUME_MAPS_', num2str(i, '%03g'))) = strcat(fName, ext);        
end
res= scLst;     oArg1= pathSollStr;



%
%
%  START:
%       [dtd, errStr]= loacal_import_SPM2_SL(fname)
%
%%%%%%%%%%%%%%%%%%%%~%%%%%%%%
function [dtd, errStr]= local_import_SPM2_SL(fName)

dtd= []; errStr= '';

pathStr= fileparts(fName);
try
    softLinkS= open(fName);
catch
    errStr= strcat(mfilename, 'local_import_SPM2_SL (error): Error during opening file');
    return
end

fieldNStr= fieldnames(softLinkS);
dataCell= struct2cell(softLinkS);
paramCell= cell(4, 3);
try
    for i= 1:3 
        %Eigenvalues
        paramCell{1, i}= fullfile(pathStr, softLinkS.(strcat('EIGVAL_', num2str(i))));
        compStr= 'XYZ';
        for ii= 1:3
            paramCell{ii + 1, i}= fullfile(pathStr, softLinkS.(strcat('EIGVECT_', num2str(i), '_', compStr(ii))));            
        end
    end
    idx= find(strncmp('VOLUME_MAPS', fieldNStr, 11));
    otherParam= cell(length(idx), 1);
    for i= 1:length(idx)
        otherParam{i}= fullfile(pathStr, softLinkS.(strcat('VOLUME_MAPS_', num2str(i, '%03g'))));
    end
catch
    errStr= strcat(mfilename, 'local_import_SPM2_SL (error): softlink file is not valid');
    return
end    

% import SPMfile
[dtd, errStr]= local_import_SPM2(paramCell(:, 1), paramCell(:, 2), paramCell(:, 3), otherParam);




%
%
%  START:
%       [eigValStrLst, errStr, eigVecStrLst, otherStrLst]= local_import_SPM2_getParam
%
%%%%%%%%%%%%%%%%%%%%~%%%%%%%%
function [e1StrLst, errStr, e2StrLst, e3StrLst, oStrLst]= local_import_SPM2_getParam

e1StrLst= {}; errStr= ''; e2StrLst= {}; e3StrLst= {}; oStrLst= {}; 

if (exist('spm_get') ~= 2) 
    errStr= strcat(mfilename, '::local_import_SPM2_getParam (error): Needed SPM02 routines weren''t found');
    return;    
end

try
    eValFiles= spm_get(3 ,'eval?_*IMAGE', {'Subject - eigenValue images'});
catch
    errStr= strcat(mfilename, '::local_import_SPM2_getParam (error): Abborted by user');
    return
end
try
    eVecFiles= spm_get(9, 'evec??_*IMAGE', {'Subject - eigenVector images'});
catch
    errStr= strcat(mfilename, '::local_import_SPM2_getParam (error): Abborted by user');
    return
end
try
    oFiles= spm_get([], '*IMAGE', {'Subject  - other images'});
catch
    oFiles= {};
end

e1StrLst= cell(4, 1); e2StrLst= cell(4, 1); e3StrLst= cell(4, 1);

%eigenVals
idx= sortcellchar(eValFiles);
e1StrLst{1}= eValFiles{idx(1)}; e2StrLst{1}= eValFiles{idx(2)}; e3StrLst{1}= eValFiles{idx(3)}; 

%eigenvectors
idx= sortcellchar(eVecFiles);
e1StrLst([2 3 4])= eVecFiles(idx([1 2 3]));
e2StrLst([2 3 4])= eVecFiles(idx([4 5 6]));
e3StrLst([2 3 4])= eVecFiles(idx([7 8 9]));

% other modalities
if ~isempty(oFiles)
    idx= sortcellchar(oFiles);
    oStrLst= oFiles(idx);
end

