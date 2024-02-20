function ea_recover_native_to_mni(patdirs)

if ischar(patdirs)
    patdirs={patdirs};
end

slicer_netstim_modules = {};
slicer_netstim_modules{end+1} = fullfile(ea_getearoot, 'ext_libs', 'SlicerNetstim', 'NetstimPreferences');

s4l = ea_slicer_for_lead;
if ~s4l.is_up_to_date()
    s4l.install();
end


for i = 1:length(patdirs)

    pt_options = ea_getptopts(patdirs{i});
    native_to_mni = [pt_options.subj.norm.transform.forwardBaseName 'ants.nii.gz'];
    mni_to_native = [pt_options.subj.norm.transform.inverseBaseName 'ants.nii.gz'];
    mni_reference = fullfile(ea_space, 't1.nii');

    if ~isfile(mni_to_native)
        warning(['Transform not available: ' mni_to_native])
    end
   
    script_source = [mfilename('fullpath') '.py'];
    script_to_run = fullfile(fileparts(pt_options.subj.norm.transform.forwardBaseName), [mfilename '.py']);
    copyfile(script_source, script_to_run)

    lines_to_append = [strcat("nativeToMNITransformFile=r'", native_to_mni, "'"),...
                       strcat("MNIToNativeTransformFile=r'", mni_to_native, "'"),...
                       strcat("MNIReferenceFile=r'", mni_reference, "'")];

    script = [char(strjoin(lines_to_append,newline)), newline, fileread(script_to_run)];
    fileID = fopen(script_to_run, 'w');
    fwrite(fileID, script, 'char');
    fclose(fileID);

    command = [' --no-splash'...
           ' --ignore-slicerrc'...
           ' --no-main-window'...
           ' --additional-module-paths "' char(strjoin(string(slicer_netstim_modules),'" "')) '"'...
           ' --modules-to-ignore DICOMDiffusionVolumePlugin'...
           ' --python-script "' script_to_run '"'];

    s4l.run(command);
    delete(script_to_run)
    if ~isfile(native_to_mni)
        warning(['Apparently something went wrong and the transformation was not recovered: ' native_to_mni]);
    else
        display(['Finished sub: ' patdirs{i}]);
    end

end
