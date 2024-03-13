function pathname = eo_uigetdir(start_path, dialog_title)
% use newer java/native window for better folder selection, with
% multiselection enbled (compared to Swing JFileChooser)
% It leverages internal java method from MATLAB(MJFileChooserPerPlatform), 
% exposing multiselection. If not available, it reverts to JFileChooser;
% based on uigetdir, uigetfile_n_dir (https://www.mathworks.com/matlabcentral/fileexchange/32555-uigetfile_n_dir-select-multiple-files-and-directories), and ea_uigetdir
% Enrico Opri, University of Michigan, 2024.

if nargin == 0 || strcmp(start_path,'') % Allow a null argument.
    start_path = pwd;
end

try
    jchooser = javaObjectEDT('com.mathworks.mwswing.MJFileChooserPerPlatform');
    is_legacy_swing=false;
catch
    jchooser = javaObjectEDT('javax.swing.JFileChooser');
    is_legacy_swing=true;
end
jchooser.setCurrentDirectory(java.io.File(start_path));

jchooser.setFileSelectionMode(javax.swing.JFileChooser.DIRECTORIES_ONLY);
if nargin > 1
    jchooser.setDialogTitle(dialog_title);
end

jchooser.setMultiSelectionEnabled(true);

if ~is_legacy_swing
    jchooser.showOpenDialog([]);
    drawnow;%force drawing to avoid Java rendering bug/hangs
    status=jchooser.getState();
else
    %javax.swing.JFileChooser has a slight different output/call
    status=jchooser.showOpenDialog([]);
    drawnow;%force drawing to avoid Java rendering bug/hangs
end
if (status == javax.swing.JFileChooser.APPROVE_OPTION)
    jFile_list=jchooser.getSelectedFiles();
    pathname{size(jFile_list, 1)}=[];
    for i=1:size(jFile_list, 1)
        pathname{i} = char(jFile_list(i).getAbsolutePath);
        if isfile(pathname{i})
            pathname{i} = fileparts(pathname{i});
        end
    end
                
elseif status == javax.swing.JFileChooser.CANCEL_OPTION
    pathname = [];
else
    error('Error occured while picking file.');
end

