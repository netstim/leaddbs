function ea_trajectory_seeg(elstruct, options)
% - Plots spheres at elstruct.coords_mm{i}
% - Builds a curved trajectory per i (uses elstruct.curve_mm{i} verbatim if present)
% - Calls ea_showelectrode(...,'plan',...) to render the electrode along that curve
%

    if ~isfield(elstruct, 'coords_mm')
        error('Input must have field coords_mm (1×N cell of [n×3]).');
    end
    
    % Sphere resolution and size
    [Xs, Ys, Zs] = sphere(12);
    r_sphere = 1.5;  % mm

    hold on
    Ncells = numel(elstruct.coords_mm);
    trajectories = cell(1, Ncells);
    plan_coords = cell(1, Ncells);
    plan_traj   = cell(1, Ncells);
    for i = 1:Ncells
        coords = elstruct.coords_mm{i};
        if isempty(coords) || size(coords,2) ~= 3
            continue
        end

        % Color for trajectory
        cmap = summer(256);
        ccol = cmap(randi(size(cmap,1)),:);
        for k = 1:size(coords,1)
            x = coords(k,1); y = coords(k,2); z = coords(k,3);
                surf(r_sphere*Xs + x, r_sphere*Ys + y, r_sphere*Zs + z, ...
                     'FaceColor', ccol, 'EdgeColor', 'none');
        end

        % Curved trajectory
        if isfield(elstruct, 'curve_mm') && numel(elstruct.curve_mm) >= i ...
                && ~isempty(elstruct.curve_mm{i}) && size(elstruct.curve_mm{i},2) == 3
            Ct = elstruct.curve_mm{i};
        elseif size(coords,1) >= 3
            P0 = coords(1,:);
            P1 = coords(end,:);
            Pm = coords(round((size(coords,1)+1)/2),:);
            A   = [0.25, 0.5; 1, 1];
            rhs = [Pm-P0; P1-P0];      
            xAB = A \ rhs;                
            a = xAB(1,:); b = xAB(2,:); c0 = P0;
            t = linspace(0,1,200).';
            Ct = a.*(t.^2) + b.*t + c0; 
            trajectories{i} = Ct;
        else
            % fallback for non-curved (straight) trajectory
            t = linspace(0,1,100).';
            Ct = (1-t).*coords(1,:) + t.*coords(end,:);
            trajectories{i} = (1-t).*coords(1,:) + t.*coords(end,:);
        end

        % Show the trajectory. Can comment out if only want to see spheres
        plot3(Ct(:,1), Ct(:,2), Ct(:,3), '-', 'LineWidth', 6, 'Color', ccol);

        % Stash for ea_showelectrode
        plan_coords{i} = coords;
        plan_traj{i}   = Ct;
    end
    hold off
    nonempty = ~cellfun(@isempty, plan_coords) & ~cellfun(@isempty, plan_traj);
    plan_coords = plan_coords(nonempty);
    plan_traj   = plan_traj(nonempty);
    if isempty(plan_coords)
        return
    end
    sides = 1:numel(plan_coords);
end