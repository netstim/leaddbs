function [elmodel,first_notempty_side]=ea_get_first_notempty_elmodel(reco_props)
    % [elmodel,first_notempty_side]=ea_get_first_notempty_elmodel(reco.props)
    % Copyright (C) 2020 Emory University, USA, School of Medicine
    % Enrico Opri

    elmodel=[];
    if ~isempty(reco_props)
        for first_notempty_side=1:length(reco_props)
            elmodel=reco_props(first_notempty_side).elmodel;
            if ~isempty(elmodel)
                break
            end
        end
    end

    %if it is still empty, return empty and throw a warning
    if isempty(elmodel)
        %no model was found
        warning('No electrode model specified.');
        elmodel=[];
    else
        % Fix Abbott lead name
        elmodel = strrep(elmodel, 'St. Jude', 'Abbott');
    end
end
