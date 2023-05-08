function ea_gencode_applynormtofile_menu(varargin)





forwardvars=varargin(3:end);
handles=varargin{3};
    uipatdir=getappdata(handles.leadfigure,'uipatdir');
forwardvars{1}=uipatdir;





[fn, pth] = uiputfile('*.m', 'Specify location for new job file...','lead_job.m');
if ~fn % user pressed cancel
    return
end

try
    options = rmfield(options, 'root');
    options = rmfield(options, 'patientname');
end

fID = fopen([pth, fn], 'w');

% export comments
fprintf(fID, '%s\n', ['function ', fn(1:end-2)]);

fprintf(fID, '%s\n',['% - Lead-DBS Job created on ', datestr(clock), ' -']);
fprintf(fID, '%s\n','% --------------------------------------');
fprintf(fID, '\n');

fprintf(fID, '%s\n',['lead path;']);
fprintf(fID, '\n');



%fprintf(fID, '%s\n',['uipatdir={'' ''}; % specify path to patient folder(s) in the following lines:']);


optionsCode=ea_gencode(forwardvars);
fprintf(fID, '%s\n', [optionsCode{1},' % Specify patient folder(s) on which to base normalizations here.']);
for e = 2:length(optionsCode)-1
    fprintf(fID, '%s\n', optionsCode{e});
end
fprintf(fID, '%s\n', [optionsCode{length(optionsCode)},' % Exchange with path to an image that defines the reference space (resolution) to be mapped to.']);

fprintf(fID, '%s\n', ['forwardvars=[forwardvars,...']);
fprintf(fID, '%s\n', ['{''/file/to/map.nii''}]; % Specify nifti file to map here. A full path to the file indicates the same file']);
fprintf(fID, '%s\n', ['%     will be used in each patient. A local path (e.g. ''anat_t2.nii'' indicates a filename inside each']);
fprintf(fID, '%s\n', ['%     patient folder will be used.']);

% add executing part:
fprintf(fID, '%s\n', ['ea_applynormtofile_menu(forwardvars{:}); % execute command.']);
fclose(fID);
edit(fullfile(pth, fn));

