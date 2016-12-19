% interface program for batch program
%
% Author: Susanne Schnell
% PC 25.07.2008

function out = dti_mapop_ui(cmd, P)

switch cmd,
    case 'normalize'
        [path1 name1] = fileparts(P.filename1{1});
        autofname = sprintf('%s_normMAP.mat',name1);
        %normalization
        probStruct = probstruct_read(P.filename1{1});
        [probStruct,errStr] = probstruct_op('NORM',probStruct);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'addition'
        [path1 name1] = fileparts(P.filename1{1});
        [path2 name2] = fileparts(P.filename2{1});
        autofname = sprintf('%s_%s_addMAP.mat',name1, name2);
        %addition
        probStruct1 = probstruct_read(P.filename1{1});
        probStruct2 = probstruct_read(P.filename2{1});
        [probStruct,errStr] = probstruct_op('ADD',probStruct1,probStruct2);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'multiplication'
        [path1 name1] = fileparts(P.filename1{1});
        [path2 name2] = fileparts(P.filename2{1});
        autofname = sprintf('%s_%s_multMAP.mat',name1, name2);
        % multiplication
        probStruct1 = probstruct_read(P.filename1{1});
        probStruct2 = probstruct_read(P.filename2{1});
        [probStruct,errStr] = probstruct_op('MULT',probStruct1,probStruct2);
        if ~isempty(errStr)
            error(errStr)
        end
    case 'connmulti'
        P1 = struct('filename1',{{}},'filename2',{{}},'newfilename','');
        out.files = {};
        for cm = 1:size(P.connmtx,1)
            P1.filename1 = P.filenames(P.connmtx(cm,1));
            P1.filename2 = P.filenames(P.connmtx(cm,2));
            out1 = dti_mapop_ui('multiplication',P1);
            out.files = [out.files(:); out1.files(:)];
        end
        return; % do not run the code below, this will be run by each
                % instance in the loop above
    case 'connmulti_check'
        if max(P.connmtx(:)) > numel(P.filenames)
            out = 'Connection matrix indices larger than number of supplied files.';
        else
            out = '';
        end
        return; % do not run the code below, this does not make sense for a
                % check function
end
% newfilename
if strcmp(P.newfilename,'.mat') || isempty(P.newfilename)
    path = fileparts(P.filename1{1});
    filename = fullfile(path, autofname);
else
    [path,name,ext] = fileparts(P.filename1{1});
    [newpath,name,ext] = fileparts(P.newfilename);
    if isempty(newpath)
        filename = fullfile(path,[name ext]);
    else
        filename = fullfile(newpath,[name ext]);
    end
end
    
%save result
ofilename = mrstruct_write(probStruct,filename);
out.files{1} = ofilename;

