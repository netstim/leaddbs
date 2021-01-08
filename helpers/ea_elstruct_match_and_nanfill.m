function elstruct=ea_elstruct_match_and_nanfill(elstruct)
    %Fill missing sides with nans, matching it to the other present side (e.g. L and R).
    %At least one side should be present, which should be always the case.
    %Enrico Opri, 2020

    %the minimum number of sides is two. First is R, second is L
    for el=1:length(elstruct)
        if isfield(elstruct(el),'coords_mm')
            if length(elstruct(el).coords_mm)==1
                elstruct(el).coords_mm{2}=[];
            elseif isempty(elstruct(el).coords_mm)
                error('need to have at least 1 tract')
            end
            %match and fill if necessary.
            elstruct(el).coords_mm=match_and_nanfill(elstruct(el).coords_mm);
        end

        if isfield(elstruct(el),'coords_acpc') && ~(~iscell(elstruct(el).coords_acpc) && isnan(elstruct(el).coords_acpc))
            if length(elstruct(el).coords_acpc)==1
                elstruct(el).coords_acpc{2}=[];
            elseif isempty(elstruct(el).coords_acpc)
                error('need to have at least 1 tract')
            end
            %match and fill if necessary.
            elstruct(el).coords_acpc=match_and_nanfill(elstruct(el).coords_acpc);
        end
        %{
        %for now do not apply this match to the trajectory
        if isfield(elstruct(el),'trajectory')
            if length(elstruct(el).trajectory)==1
                elstruct(el).trajectory{2}=[];
            elseif isempty(elstruct(el).trajectory)
                error('need to have at least 1 tract')
            end
            %match and fill if necessary.
            elstruct(el).trajectory=match_and_nanfill(elstruct(el).trajectory);
        end
        %}   
    end
end

function coords=match_and_nanfill(coords)
    num_points=nan;%contacts in the case of coords_mm and coords_acpc, points in the case of trajectory
    num_axes=nan;
    for iside=1:length(coords)
        side=iside;%just for readability, as this is the side
        if ~isempty(coords{side})
            num_points=size(coords{side},1);
            num_axes=size(coords{side},2);
            break;
        end
    end
    if isnan(num_points) || isnan(num_axes)
        error('This should not happen, is there no electrode set yet? Remember to run the electrode/lead reconstruction first (Panel 5)');
    end

    for iside=1:length(coords)
        %side=options.sides(iside);
        side=iside;%just for readability, as this is the side
        if isempty(coords{side})
            %fill with nans, to leave a placeholder for no contact/electrode/lead
            coords{side}=nan(num_points,num_axes);
        end
    end
end