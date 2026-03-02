function plot2datlasstruct(path, atlas, plane, coord, thr) 

    for side = 1:2
        for roi = 1: length(atlas.names)
            if side ==1
                 nii = ea_load_nii([path filesep 'rh' filesep atlas.names{roi}]); 
            else
                 nii = ea_load_nii([path filesep 'lh' filesep atlas.names{roi}]); 
            end

            switch plane
                case 'x'
                    x_coord = coord; 
                    
                    % find idx of coord
                    x_idx = round([x_coord, 1, 1,1]/nii.mat'); 
                    x_idx = x_idx(1);
                    
                    if x_idx>1 && x_idx<size(nii.img, 1) % skip if the structure does not cross the slice
                        slice = nii.img(x_idx, :, :);
                        
                        y_coords =  [ones(1,size(slice, 2)); 1:size(slice, 2); ones(1,size(slice, 2)); ...
                            ones(1,size(slice, 2))]' *nii.mat';
                        y_coords = y_coords(:, 2); 
                        z_coords =  [ones(1,size(slice, 3)); ones(1,size(slice, 3)); 1:size(slice, 3); ...
                            ones(1,size(slice, 3))]' *nii.mat';
                        z_coords = z_coords(:, 3); 
                        
                        if 0 % plots to check things
                            figure
                            subplot(121)
                            imagesc(y_coords, z_coords', squeeze(slice))
                            xlabel('y'), ylabel('z'), axis square, axis xy
                        end 
          
                        y_coords_upsampled = linspace(y_coords(1), y_coords(end), size(y_coords, 1)*2); 
                        z_coords_upsampled = linspace(z_coords(1), z_coords(end), size(z_coords, 1)*2); 
                        slice_upsampled = interp2(y_coords', z_coords', squeeze(slice)', ...
                            y_coords_upsampled', z_coords_upsampled, 'cubic');
                        
                        if 0 % plots to check things
                            subplot(1,2,2)
                            imagesc(y_coords_upsampled, z_coords_upsampled, slice_upsampled)
                            xlabel('y'), ylabel('z'), axis square, axis xy
                        end 
        
                        c = contourc(y_coords_upsampled, z_coords_upsampled, slice_upsampled, [thr thr]);
                        c = c(:, 2:end); 
                        if 0
                            hold on
                            plot(c(1,:), c(2,:), 'linew', 2, 'color', [1 0 0])
                        end
                        
                        % ea_mnifigure
                        plot3(ones(1,length(c)-1)*x_coord, c(1,2:end), c(2,2:end), 'linew', 2, 'color', atlas.colormap(roi, :))
                        x_patch = ones(1, size(c,2)-1) * x_coord;
                        y_patch = c(1, 2:end);
                        z_patch = c(2, 2:end);

                        patch(x_patch, y_patch, z_patch, atlas.colormap(roi, :), ...
                            'FaceAlpha', 0.3, ...         % adjust transparency (0 = fully transparent, 1 = opaque)
                            'EdgeColor', 'none');         % no edge line

                    end
                
                case 'y' 
                    y_coord = coord; 
                    
                    % find idx of coord
                    y_idx = round([1, y_coord, 1,1]/nii.mat'); 
                    y_idx = y_idx(2);
                    
                    if y_idx>1 && y_idx<size(nii.img, 2) % skip if the structure does not cross the slice
                        slice = nii.img(:, y_idx, :);
                        
                        x_coords =  [1:size(slice, 1); ones(1,size(slice, 1)); ones(1,size(slice, 1)); ...
                            ones(1,size(slice, 1))]' *nii.mat';
                        x_coords = x_coords(:, 1); 
                        z_coords =  [ones(1,size(slice, 3)); ones(1,size(slice, 3)); 1:size(slice, 3); ...
                            ones(1,size(slice, 3))]' *nii.mat';
                        z_coords = z_coords(:, 3); 
                        
                        if 0 % plots to check things
                            figure
                            subplot(121)
                            imagesc(x_coords, z_coords', squeeze(slice))
                            xlabel('x'), ylabel('z'), axis square, axis xy
                        end 
          
                        x_coords_upsampled = linspace(x_coords(1), x_coords(end), size(x_coords, 1)*2); 
                        z_coords_upsampled = linspace(z_coords(1), z_coords(end), size(z_coords, 1)*2); 
                        slice_upsampled = interp2(x_coords', z_coords', squeeze(slice)', ...
                            x_coords_upsampled', z_coords_upsampled, 'cubic');
                        
                        if 0 % plots to check things
                            subplot(1,2,2)
                            imagesc(x_coords_upsampled, z_coords_upsampled, slice_upsampled)
                            xlabel('x'), ylabel('z'), axis square, axis xy
                        end 
        
                        c = contourc(x_coords_upsampled, z_coords_upsampled, slice_upsampled, [thr thr]);
                        c = c(:, 2:end); 
                        if 0
                            hold on
                            plot(c(1,:), c(2,:), 'linew', 2, 'color', [1 0 0])
                        end
                        
                        % ea_mnifigure
                        plot3(c(1,2:end), ones(1,length(c)-1)*y_coord,  c(2,2:end), 'linew', 2, 'color', atlas.colormap(roi, :))
                        x_patch = c(1, 2:end);
                        y_patch = ones(1, size(c,2)-1) * y_coord;
                        z_patch = c(2, 2:end);

                        patch(x_patch, y_patch, z_patch, atlas.colormap(roi, :), ...
                            'FaceAlpha', 0.3, ...         % adjust transparency (0 = fully transparent, 1 = opaque)
                            'EdgeColor', 'none');         % no edge line

                    end                

                case 'z'
                    z_coord = coord; 
                    
                    % find idx of coord
                    z_idx = round([1,1,z_coord,1]/nii.mat'); 
                    z_idx = z_idx(3);
                    
                    if z_idx>1 && z_idx<size(nii.img, 3) % skip if the structure does not cross the slice
                        slice = nii.img(:, :, z_idx);
                        
                        x_coords =  [1:size(slice, 1); ones(1,size(slice, 1)); ones(1,size(slice, 1)); ...
                            ones(1,size(slice, 1))]' *nii.mat';
                        x_coords = x_coords(:, 1); 
                        y_coords =  [ones(1,size(slice, 2)); 1:size(slice, 2); ones(1,size(slice, 2)); ...
                            ones(1,size(slice, 2))]' *nii.mat';
                        y_coords = y_coords(:, 2); 
                        
                        if 0 % plots to check things
                            figure
                            subplot(121)
                            imagesc(x_coords, z_coords', squeeze(slice))
                            xlabel('x'), ylabel('z'), axis square, axis xy
                        end 
          
                        x_coords_upsampled = linspace(x_coords(1), x_coords(end), size(x_coords, 1)*2); 
                        y_coords_upsampled = linspace(y_coords(1), y_coords(end), size(y_coords, 1)*2); 
                        slice_upsampled = interp2(x_coords', y_coords', squeeze(slice)', ...
                            x_coords_upsampled', y_coords_upsampled, 'cubic');
                        
                        if 0 % plots to check things
                            subplot(1,2,2)
                            imagesc(x_coords_upsampled, z_coords_upsampled, slice_upsampled)
                            xlabel('x'), ylabel('z'), axis square, axis xy
                        end 
        
                        c = contourc(x_coords_upsampled, y_coords_upsampled, slice_upsampled, [thr thr]);
                        c = c(:, 2:end); 
                        if 0
                            hold on
                            plot(c(1,:), c(2,:), 'linew', 2, 'color', [1 0 0])
                        end
                        
                        % ea_mnifigure
                        plot3(c(1,2:end),c(2,2:end),ones(1,length(c)-1)*z_coord,   'linew', 2, 'color', atlas.colormap(roi, :))
                        x_patch = c(1, 2:end);
                        y_patch = c(2, 2:end);
                        z_patch = ones(1, size(c,2)-1) * z_coord;

                        patch(x_patch, y_patch, z_patch, atlas.colormap(roi, :), ...
                            'FaceAlpha', 0.3, ...         % adjust transparency (0 = fully transparent, 1 = opaque)
                            'EdgeColor', 'none');         % no edge line

                    end     
            end % switch plane
        end % for roi 
    end % for side
end % function