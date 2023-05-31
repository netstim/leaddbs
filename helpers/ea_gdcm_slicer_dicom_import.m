function ea_gdcm_slicer_dicom_import(dicom_dir, output_dir)

script_source = [mfilename('fullpath') '.py'];
script_to_run = fullfile(output_dir, [mfilename '.py']);
copyfile(script_source, script_to_run)

lines_to_append = [strcat("dicomDataDir='", dicom_dir, "'"),...
                   strcat("outputDir='", output_dir, "'")];

script_read = fileread(script_to_run);
script_read = [char(strjoin(lines_to_append,newline)), newline, script_read];
fileID = fopen(script_to_run, 'w');
fwrite(fileID, script_read, 'char');
fclose(fileID);

command = [' --no-splash'...
       ' --ignore-slicerrc'...
       ' --no-main-window'...
       ' --python-script "' script_to_run '"'];

s4l = ea_slicer_for_lead;
if ~s4l.is_installed_and_up_to_date()
    s4l.install();
end

s4l.run(command)

end

