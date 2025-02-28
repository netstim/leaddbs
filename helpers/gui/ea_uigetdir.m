function pathname = ea_uigetdir(start_path, dialog_title);
    if ~exist('start_path','var')
        start_path=pwd;
    end
    if exist('dialog_title','var')
        try
            pathname = eo_uigetdir(start_path, dialog_title);
        catch
        end
    else
        pathname = eo_uigetdir(start_path);
    end
end