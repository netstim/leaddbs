% interface program for batch program
% dti_read_data_ui: user interface for reading the raw data (main: dti_cfg_read.m)
%
% Author: Susanne Schnell
% PC 10.04.2008

function out = dti_readdata_ui(cmd,P)

switch cmd
    case 'binary'
        [path,nam,ext] = fileparts(P.binary{1});
        if strcmp(ext,'.bin')
            % open data file and check if all necessary diffusion
            % information provided
            [res, err] = dw_data_admin('open',P.binary{1}(1:end-8));
            if isempty(res)
                error(err)
            end
            if res == 0
                error(err)
            end
            DE_scheme = dw_data_admin('get_diffDirs');
            bvalue = dw_data_admin('get_bVal');
            nob0s = dw_data_admin('get_user');
            dw_data_admin('close');
            if ~isempty(nob0s)
                nob0s = nob0s.nob0s;
            end
            if ~isempty(DE_scheme) && ~isempty(bvalue) && ~isempty(nob0s)
                out.files{1} = fullfile(path, [nam,ext]);
            else
                out.files = {};
                error('The necessary diffusion information is missing. This data cannot processed using the batch modus. Please use the DTI&FiberTools GUI instead.');
            end
        else
            out.files = {};
            error('This input data file type is not supported yet.');
        end
    case 'mrstruct'
        [path,nam,ext] = fileparts(P.filename{1});
        if strcmp(ext,'.mat')
            mr = mrstruct_read(P.filename{1});
            mr.user.nob0s = P.nob0s;
            mr.user.bvalue = P.bvalue;
            mr.user.DEscheme = P.descheme;
            mrstruct_write(mr,P.filename{1});
            out.files{1} = P.filename{1};
        else
            out.files = {};
            error('This input data file type is not a matlab variable and is not supported yet.');
        end
    case 'append'
        [path,nam,ext] = fileparts(P.filename2{1});
        if strcmp(ext,'.mat')
            mr = mrstruct_read(P.filename2{1});
            dtd = dtdstruct_read(P.dtdname{1});
            [dtd, err] = dtdstruct_modify(dtd,'addStruct',mr,P.volname);
            if ~isempty(err)
                error(err);
                out.files = {};
            else
                dtdstruct_write(dtd,P.dtdname{1});
                out.files{1} = P.dtdname{1};
            end
        else
            out.files = {};
            error('This input data file type is not a matlab variable and is not supported yet.');
        end
    case 'dicom'
        numberofdirs = numel(P.dicom);
        for m = 1 : numberofdirs
            [path,nam,ext] = fileparts(P.dicom{m});
            list = dir(fullfile(path,'*.dcm'));
            if m == 1
                [out.files{1}, status, errorstring] = read_dti_dicom(path, list,0,[],1);
            else
                [out.files{1}, status, errorstring] = read_dti_dicom(path, list,1,out.files{1},m);
            end
            if ~isempty(errorstring)
                error(errorstring);
            end
            if ~isequal(status,[1 1 1])
                fprintf('The Dicom data does not privide all necessary information for the calculation of the tensors.\nThe reason for this might be that you don not work with Siemens Dicom.\n Please use the GUI instead (dti_tool) of the Batch Editor (see manual).\nAnother option would be that you read in your data into an mrstruct (see manual) and set up the batch using this as raw data.');
                return
            end
        end
    case 'bruker'
        [out.files{1}, status, errorstring] = bruker_read_DTI_data('Filepath',P.bruker{1});
        if ~isempty(errorstring)
            error(errorstring);
        end
        if ~isequal(status,[1 1 1])
            error('The Bruker data does not privide all necessary information for the calculation of the tensors.');
        end       
    case 'computeDTDfromHARDI'        
        hr = mrstruct_read(P.filename{1});
        
        if isfield(P.newfile,'auto')
            [p n e] = fileparts(P.filename{1});
            if not(isempty(strfind(n,'HARDI'))),
                n = strrep(n,'HARDI','DTD');
            else
                n = [n '_DTD'];
            end
            outname = fullfile(p,[n e]);
        else
            outname = fullfile(P.newfile.out.dir{1}, P.newfile.out.fnameDTD);
        end;
        out.files{1} = outname;     
        kurto = false;
        if isfield(P.kurto,'kyes'),
            kurto = true;
        end;
        dtd = convertHARDI2DTD(hr,P.threshold,kurto);
        [res, errStr] = dtdstruct_write(dtd,outname);        
        if ~isempty(errStr)
            error(errStr)
        end        
        return;
        
end
