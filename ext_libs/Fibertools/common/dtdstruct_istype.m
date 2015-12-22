function [res, typeNameStr, errStr]= dtdstruct_istype(dtdStruct)
%function [res, typeNameStr, errStr]dtdstruct_istype(arg)
%
%   function determines if arg is a valid dtdStruct. If this is the case,
%   also the type name is returned as the second return value. 
%
%   arg:    arbitary data .. may be a dtdStruct
%
%   return values
%   res:    1 if arg is an dtdStrcut 0 else
%   typeNameStr: String identifying the dtdStruct type
%   errStr: If an error occured, errStr is identify the error
%
% Bjoern W. Kreher
% 08/02
%
% UNIX


res= 0; errStr= ''; typeNameStr= '';
if isstruct(dtdStruct)
    
    if mrstruct_istype(dtdStruct)
        if length(mrstruct_query(dtdStruct, 'sizeAy')) == 3
            res= 1;
            typeNameStr= 'mrStruct';
        else
            errStr= sprintf('%s(warning): size of mrStructhas wrong dimension for dtdStruct conversion', mfilename);
        end
        return;
    end
    
    dataCell= struct2cell(dtdStruct);
    fieldStr= fieldnames(dtdStruct);
    if size(dataCell) == 1
        return
    end
    
    %%% test ob soft link struct
    ok= 1;
    for i= 1:length(dataCell)
        ok= ok & ischar(dataCell{i});
    end
    if ok && ~isempty(find(strcmp(fieldStr, 'PRIV_FORMAT'), 1))
        res= 1;     typeNameStr= 'soft-link dtdStruct';
        return
    end
    
    ok= 1;
    sizeSoll= [];
    dimsAy= zeros(length(dataCell), 1);
    for i= 1:length(dataCell)
        if mrstruct_istype(dataCell{i}) 
            sizeAy= mrstruct_query(dataCell{i}, 'sizeAy');
            if isempty(sizeSoll)
                if length(sizeAy) < 3
                    sizeSoll= ones(1, 3);
                    sizeSoll(1:length(sizeAy))= sizeAy;
                else
                    sizeSoll= sizeAy(1:3);
                end
            end
%             if (length(sizeAy) < 3) || ~isequal(sizeAy(1:3), sizeSoll(1:3))
%                 ok= 0;
%                 errStr= 'dtdstruct_istype (error): Data modalities of different size';
%             end
            dimsAy(i)= length(sizeAy); %mrstruct_query(dataCell{i}, 'dimensions');
        else
            if ~strcmp(fieldStr{i}, 'version')
                ok= 0;
            end
            dimsAy(i)= 3;
        end
    end
    if (ok == 1) && (isempty(find(dimsAy ~= 3, 1)))
        typeNameStr= 'someMR';
        res= 1;
        return
    elseif isfield(dtdStruct, 'velocityVect_struc')
        res= 1;
        typeNameStr= 'VelVD';
        
    elseif isfield(dtdStruct, 'eigenVal_struc') && isfield(dtdStruct, 'eigenVec_struc')
        res= 1;
        typeNameStr= 'DTD';
        
        if isfield(dtdStruct, 'eigVal_1') && isfield(dtdStruct, 'eigVec_1') && isfield(dtdStruct, 'eigVal_3') && isfield(dtdStruct, 'p_Value')
            typeNameStr= 'monoMDT';
            res= 1;
            if isfield(dtdStruct, 'eigVal_2') && isfield(dtdStruct, 'eigVec_2')
                typeNameStr= 'multiDTD';
                res= 1;
            end
        end
    end
else
    errStr= 'dtdstruct_istype (error): Argument is not a struct';
end

if res == 0
    errStr= 'dtdstruct_istype (error): data ist not a dtdStruct';
end