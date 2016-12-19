% interface program for batch modus
%
% Author: Susanne Schnell
% PC 10.04.2008

function out = dti_tensor_ui(P)

if nargin == 1 && isstruct(P) && isfield(P, 'filename')
    [path,nam,ext] = fileparts(P.filename{1});
    filename = fullfile(path,nam);
    savemode = bitset(bitset(0,1,P.singleslice),2,P.hardi);
    if strcmp(ext,'.bin')
        [ok, msg] = calculate_dti(filename(1:end-4),P.threshold,savemode);
        out.files{1} = strcat(filename(1:end-4),'_DTD.mat');
    elseif strcmp(ext,'.mat')
        mr = mrstruct_read(strcat(filename,'.mat'));
        if isfield(mr.user,'bvalue') && isfield(mr.user,'DEscheme') && isfield(mr.user,'nob0s')
            if iscell(mr.user.DEscheme)
                if exist(mr.user.DEscheme{1},'file')
                    [ok, msg] = calculate_dti(strcat(filename,'.mat'),P.threshold,savemode,load(mr.user.DEscheme{1}),mr.user.bvalue,mr.user.nob0s);
                    out.files{1} = strcat(filename,'_DTD.mat');
                else
                    [ok, msg] = calculate_dti(strcat(filename,'.mat'),P.threshold,savemode,mr.user.DEscheme,mr.user.bvalue,mr.user.nob0s);
                    out.files{1} = strcat(filename,'_DTD.mat');
                end
            else
                if exist(mr.user.DEscheme,'file')
                    [ok, msg] = calculate_dti(strcat(filename,'.mat'),P.threshold,savemode,load(mr.user.DEscheme),mr.user.bvalue,mr.user.nob0s);
                    out.files{1} = strcat(filename,'_DTD.mat');
                else
                    [ok, msg] = calculate_dti(strcat(filename,'.mat'),P.threshold,savemode,mr.user.DEscheme,mr.user.bvalue,mr.user.nob0s);
                    out.files{1} = strcat(filename,'_DTD.mat');
                end
            end
        else
            error('Mrstruct does not contain the necessary information in field user (bvalue, DEscheme or nob0s is missing)');
        end
        clear mr
    elseif isempty(ext)
        if strcmp(filename(end-3:end),'_raw')
            [ok, msg] = calculate_dti(filename(1:end-4),P.threshold,savemode);
            out.files{1} = strcat(filename(1:end-4),'_DTD.mat');
        else
            [ok, msg] = calculate_dti(filename,P.threshold,savemode);
            out.files{1} = strcat(filename,'_DTD.mat');
        end
    else
        error('The filename for the raw data file is wrong or not supported.');
    end
    if ~strcmp(msg,'DTI calculation finished.')
        error(msg);
    end
end

