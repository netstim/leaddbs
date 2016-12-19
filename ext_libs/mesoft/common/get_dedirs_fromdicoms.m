
function [dedirs,bval] = get_dedirs_fromdicoms

path = uigetdir(pwd,'Select folder with dicoms');

if ispc
    files = dir(strcat(path,'\*.dcm'));
else
    files = dir(strcat(path,'/*.dcm'));
end

% extract diffusion sequence specific information, seems to work for
% Siemens VA25 and VB12 and VB13!
% edges müssen existieren, um die Diffusionsrichtungen in den Patientenraum
% zu drehen
for m = 1 : size(files,1)
    dirs = []; b_value = [];
    header(m) = dicom_read_singlefile(fullfile(path,files(m,:).name),1,[],1);
    if isfield(header(m),'Private_0029_1010')
        for n = 1 : length(header(m).Private_0029_1010)
            if strcmp(header(m).Private_0029_1010(n).name,'DiffusionGradientDirection')
                if str2num(header(m).Private_0029_1010(n).item(1,1).val) ~= 0
                    dirs = [str2num(header(m).Private_0029_1010(n).item(1,1).val),str2num(header(m).Private_0029_1010(n).item(1,2).val),str2num(header(m).Private_0029_1010(n).item(1,3).val)];
                else
                    dirs = [0 0 0];
                end
            end
            if strcmp(header(m).Private_0029_1010(n).name,'B_value' )
                b_value = str2num(header(m).Private_0029_1010(n).item(1,1).val);
            end
            if ~isempty(b_value) && ~isempty(dirs)
                bval(m) = b_value;
                dedirs(m,:,:,:) = dirs;
                break;
            end
        end
    end
end
bval = bval';