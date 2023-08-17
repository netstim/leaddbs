function [min_bound_per_contact, max_bound_per_contact, S] = ea_get_currents_per_contact(min_conc,max_conc,min_segm,max_segm, reconst, side, check_S)

    % create bounds (limits) for current optimization
    % and adjust the amplitude in S

    % min_conc,max_conc,min_segm,max_segm - thresholds for contacts
    % reconst - lead-dbs electrode reconstruction file
    % side - (OSS-DBS notation, 0 = right_hemisphere, 1 = left_hemisphere)
    % check_S - if 1, adjust amplitude on sources

    if size({reconst.reco.props.elmodel},2) == 2
        %sides = 2;  % both sides implanted
        el_models{1} = strrep(reconst.reco.props(1).elmodel, ' ','_');
        el_models{2} = strrep(reconst.reco.props(2).elmodel, ' ','_');
        if ~strcmp(el_models{1},el_models{2})
            disp("Bilateral stim with different electrode models currently not supported")
            return
        end
    elseif reconst.reco.mni.markers(1).head(1) > 0.0
        %sides = 0; % right side
        el_models{1} = strrep(reconst.reco.props(1).elmodel, ' ','_');
        el_models{2} = -1;
    else
        %sides = 1; % left side
        el_models{2} = -1;
        el_models{1} = strrep(reconst.reco.props(1).elmodel, ' ','_');
    end

    % electrode models, should be expanded
    concentric4 = ['PINS_Medical_L303','PINS_Medical_L302','PINS_Medical_L301', 'Medtronic_3389','Medtronic_3387','Medtronic_3391','St._Jude ActiveTip_(6142-6145)', 'St._Jude_ActiveTip_(6146-6149)','Abbott ActiveTip (6146-6149)','Abbott ActiveTip (6142-6145)'];
    concentric8 = ['Boston_Scientific_Vercise'];
    segmented8 = ['St._Jude_Directed_6172_(short)', 'St._Jude_Directed_6180','St._Jude_Directed_6173_(long)','Boston_Scientific_Vercise_Directed','Abbott Directed 6172 (short)','Abbott Directed 6173 (long)', 'Medtronic B33005', 'Medtronic B33015'];
    
    % add a sanity check for current limits

    % use current thresholds to set amplitudes 
    if any(contains(concentric4,el_models{1}))
        min_bound_per_contact = ones(1,4) * min_conc;
        max_bound_per_contact = ones(1,4) * max_conc;

        if check_S == 1
            % we need to change S for 4-electrode
            load([ea_getearoot,'cleartune/cleartunePAM/S_4contact.mat']);
            if min_conc < 0.0
                if length(side) == 2
                    S.Rs1.amp = 4 * (max_conc - min_conc);
                    S.Ls1.amp = 4 * (max_conc - min_conc);
                elseif side == 0
                    S.Rs1.amp = 4 * (max_conc - min_conc); 
                else
                    S.Ls1.amp = 4 * (max_conc - min_conc);
                end
            else
                if length(side) == 2
                    S.Rs1.amp = 4 * (max_conc);
                    S.Ls1.amp = 4 * (max_conc);
                elseif side == 0
                    S.Rs1.amp = 4 * (max_conc);
                else
                    S.Ls1.amp = 4 * (max_conc);
                end   
            end
        end
    elseif any(contains(concentric8,el_models{1}))
        min_bound_per_contact = ones(1,8) * min_conc;
        max_bound_per_contact = ones(1,8) * max_conc;

        if check_S == 1
            % we need to change S for 8-electrode same %
            load([ea_getearoot,'cleartune/cleartunePAM/S_8contact_conc.mat']);
            if min_conc < 0.0
                if length(side) == 2
                    S.Rs1.amp = 8 * (max_conc - min_conc);
                    S.Ls1.amp = 8 * (max_conc - min_conc);
                elseif side == 0
                    S.Rs1.amp = 8 * (max_conc - min_conc);
                else
                    S.Ls1.amp = 8 * (max_conc - min_conc);
                end
            else
                if length(side) == 2
                    S.Rs1.amp = 8 * (max_conc);
                    S.Ls1.amp = 8 * (max_conc);
                elseif side == 0
                    S.Rs1.amp = 8 * (max_conc);
                else
                    S.Ls1.amp = 8 * (max_conc);
                end   
            end
        end
    elseif any(contains(segmented8,el_models{1}))
        min_bound_per_contact = ones(1,8);
        min_bound_per_contact(1,1) = min_conc;
        min_bound_per_contact(1,8) = min_conc;
        min_bound_per_contact(1,2:7) = min_segm;

        max_bound_per_contact = ones(1,8);
        max_bound_per_contact(1,1) = max_conc;
        max_bound_per_contact(1,8) = max_conc;
        max_bound_per_contact(1,2:7) = max_segm;

        if check_S == 1
            load([ea_getearoot,'cleartune/cleartunePAM/S_8contact_segm.mat']);
            if min_conc < 0.0  % we assume that both concentric and segmented will have either both or one polarity
                if length(side) == 2
                    S.Rs1.amp = 2 * (max_conc - min_conc) + 6 * (max_segm - min_segm);
                    S.Ls1.amp = 2 * (max_conc - min_conc) + 6 * (max_segm - min_segm);
                elseif side == 0
                    S.Rs1.amp = 2 * (max_conc - min_conc) + 6 * (max_segm - min_segm);
                else
                    S.Ls1.amp = 2 * (max_conc - min_conc) + 6 * (max_segm - min_segm);
                end
            else
                if length(side) == 2
                    S.Rs1.amp = 2 * (max_conc) + 6 * (max_segm);
                    S.Ls1.amp = 2 * (max_conc) + 6 * (max_segm);
                elseif side == 0
                    S.Rs1.amp = 2 * (max_conc) + 6 * (max_segm);
                else
                    S.Ls1.amp = 2 * (max_conc) + 6 * (max_segm);
                end   
            end
        end
    else
        ea_warndlg('The electrode model is not recongnized')
        return
    end

    if check_S == 1
        if length(side) == 2
            S.amplitude{1,1}(1) = S.Rs1.amp;
            S.amplitude{1,2}(1) = S.Ls1.amp;
        elseif side == 0
            S.amplitude{1,1}(1) = S.Rs1.amp;
            S.amplitude{1,2}(1) = 0.0; % one side at a time
        else
            S.amplitude{1,2}(1) = S.Ls1.amp;
            S.amplitude{1,1}(1) = 0.0; % one side at a time
        end
    else
        S = 0;  % won't be used
    end
end