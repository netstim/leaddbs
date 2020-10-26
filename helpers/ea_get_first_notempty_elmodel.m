function elmodel=ea_get_first_notempty_elmodel(reco___props)
    % elmodel=ea_get_first_notempty_elmodel(reco.props)
    % Copyright (C) 2020 Emory University, USA, School of Medicine
    % Enrico Opri

    elmodel=[];
    if ~isempty(reco___props)
        for side=1:length(reco___props)
            elmodel=reco___props(side).elmodel;
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
    end
end