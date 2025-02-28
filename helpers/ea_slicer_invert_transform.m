function [] = ea_slicer_invert_transform(input_transform, reference_volume, output_transform)

% Set-up Custom Slicer
s4l = ea_slicer_for_lead;
if ~s4l.is_up_to_date()
  s4l.install();
end

python_script = [mfilename("fullpath") '.py'];
slicer_cmd = {'--no-splash', '--no-main-window', '--ignore-slicerrc', '--python-script', ...
              ea_path_helper(python_script), ...
              ea_path_helper(input_transform), ...
              ea_path_helper(reference_volume), ...
              ea_path_helper(output_transform)};
status = s4l.run(strjoin(slicer_cmd, ' '));
if status ~= 0
    ea_error('Failed to invert the transformation!', showdlg=false, simpleStack=true);
end

end