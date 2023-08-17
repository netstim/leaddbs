function [reconst, el_model_right, el_model_left, sides] = ea_get_reconstruction(reconstruction_file_path)

    try
        reconst = load(reconstruction_file_path);
    catch
        ea_warndlg('The reconstruction file was not found')
        return
    end
    
    if length(reconst.reco.props) == 2
        sides = [0,1];  % both sides implanted
        el_model_right = strrep(reconst.reco.props(1).elmodel, ' ','_');
        el_model_left = strrep(reconst.reco.props(2).elmodel, ' ','_');
    elseif reconst.reco.mni.markers(1).head(1) > 0.0
        sides = 0; % right side
        el_model_right = strrep(reconst.reco.props(1).elmodel, ' ','_');
        el_model_left = 'None';
    else
        sides = 1; % left side
        el_model_right = 'None';
        el_model_left = strrep(reconst.reco.props(1).elmodel, ' ','_');
    end
end
