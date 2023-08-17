function [coords_mm,trajectory,markers] = ea_refinecoords(options)
%% Refine the fiducial markers following TRAC/CORE reconstruction
%  Last Revision: 18/05/2018
%  Thushara Perera (c) 2018 Bionics Institute
%  Input:
%   - lead dbs options struct
%  Output:
%   - refined reconstruction is saved to file using ea_save_reconstruction()
%   - figures useful for de-bugging displayed when is_debug = 1
%
%  Head and tail fiducial markers are repositioned on to local maxima
%  (hyperintese portion of the CT). At present refinement using post-op MRI
%  has not been implemented. Only working for Medtronic 3387 and 3389
%  electrodes.
%
%  Bits and pieces of code have been copyied from other areas of Lead-DBS.
%  There might be a better way to implement what I'm trying to do and
%  several things are yet to be implemented.
    disp('Refining fiducials by shifting them to local maxima.');

    is_debug = 0;
    saveimg = 0;
    can_export = 0;
    sample_width = 10;
    doxx = 1; % TODO: Try it in the other planes as well??

    options.prefs = ea_prefs;
    if isfield(options.prefs.reco, 'saveimg')
        saveimg = options.prefs.reco.saveimg;
    end
    if isfield(options.prefs.reco, 'exportfiducials')
        can_export = ischar(options.prefs.reco.exportfiducials);
    end
    if ~saveimg
        disp('Set prefs.reco.saveimg=1 in preference file to save marker visualisation as image.');
    end
    if ~can_export
        disp('Set prefs.reco.exportfiducials in preference file to export fiducial markers as CSV file.');
    end

    options.native = 1;
    options.loadnativereco = 1; % Load native reco intead of scrf
    [coords_mm,trajectory,markers,elmodel,~]=ea_load_reconstruction(options);

    switch options.subj.postopModality
        case 'MRI'
            disp('ea_refinecoords not implemented yet for MRI. Skipping...');
            return;
        case 'CT'
            V = spm_vol(options.subj.coreg.anat.postop.CT);
    end


    for side = options.sides
        options.elside = side;
        meantrajectory = genhd_inside(trajectory{side});
        imat = ea_resample_planes(V, meantrajectory', sample_width, doxx, 0.2);

        if doxx % span surface in x direction
            spanvector=[sample_width,0,0];
        else % span surface in y direction
            spanvector=[0,sample_width,0];
        end

        boundingbox=[meantrajectory(1,:)-spanvector;...
            meantrajectory(1,:)+spanvector;...
            meantrajectory(end,:)-spanvector;...
            meantrajectory(end,:)+spanvector];

        xx=[boundingbox(1,1),boundingbox(2,1);boundingbox(3,1),boundingbox(4,1)];
        yy=[boundingbox(1,2),boundingbox(2,2);boundingbox(3,2),boundingbox(4,2)];
        zz=[boundingbox(1,3),boundingbox(2,3);boundingbox(3,3),boundingbox(4,3)];

        % threshold slice to remove background
        b = imat;
        b(b<(max(imat(:))*0.2)) = 0; % 0.2 arbitrarily chosen based on trial and error

        % find voxel dimensions in mm
        deltax = (xx(1,2)-xx(1,1))/size(b,2);
        deltay = (yy(2,1)-yy(1,1))/size(b,1);
        deltaz = (zz(2,1)-zz(1,1))/size(b,1);
        delta = (sqrt((yy(2,1)-yy(1,1))^2 + (zz(2,1)-zz(1,1))^2))/size(b,1); % voxel dimension along trajectory

        % define inter-electrode distance
        switch options.elmodel % TODO: Set for other electrode models
            case 'Medtronic 3387'
                elgap = floor(3/delta); % 3mm inter-electrode distance
            case 'Medtronic 3389'
                elgap = floor(2/delta);
            otherwise
                disp(['Inter-electrode distance not set for ', options.elmodel, '. Skipping...']);
                return;
        end

        % filter and find first peak (head fiducial marker)
        filtered_max = ea_sgolayfilt(max(b, [], 2), 1, 21);
        if (deltaz < 0) % sometimes the slice is upside down
            filtered_max = flipud(filtered_max);
        end

        [~, idy] = ea_findpeaks(filtered_max, 'NPeaks', 1);
        if isempty(idy)
            warning(['Could not find head of electrode. Check trajectory. Skipping side ', num2str(side)]);
            continue;
        end
        if (deltaz < 0)
            idy = size(b,1) - idy;
        end
        pidy = [idy; idy-elgap*3]; % set the tail based on electrode height (FIX: for electrodes with more than 4 contacts)

        % define lateral location of fiducial
        [~, idxx] = max(imat, [], 2); % better results obtained when using orignial slice (without threshold)
        idxx = round(mean(idxx));
        pidx = [idxx; idxx];
        % [~, pidx] = max(imat(pidy,:), [], 2); % didn't work very well

        % find centre of artefact and add offset accordingly
        [cx, cy] = find_centre(b, pidx, pidy, ceil(abs(1.5/deltax)), is_debug); % 1.5mm electrode height for Medtronic
        cx(isnan(cx)) = 0;
        cy(isnan(cy)) = 0;
        pidx(1) = pidx(1) + cx(1); % only moving head laterally, difficult to locate centre of tail consistently
        pidy = pidy + cy;

        % check for errors
        if any(isnan(pidx)) || any(isnan(pidy))
            warning(['Could not find one or more fiducial. Check trajectory. Skipping side ', num2str(side)]);
            continue;
        end

        % convert to mm
        [t_head, t_tail] = convert2mm(pidx, pidy, deltax, deltay, deltaz, xx, yy, zz, b);

        % refine fiducials using axial slices
        [ahx, ahy, head_slice, hx, hy] = refine_axial(V, t_head, is_debug);
        [atx, aty, tail_slice, tx, ty] = refine_axial(V, t_tail, is_debug);
        pidx(1) = pidx(1) + ahx;
        pidx(2) = pidx(2) + atx;
        pidy(1) = pidy(1) + ahy;
        pidy(2) = pidy(2) + aty;

        % convert to mm
        [head, tail] = convert2mm(pidx, pidy, deltax, deltay, deltaz, xx, yy, zz, b);

        % apply same offset as head to x and y
        markers(side).x = markers(side).x + head - markers(side).head;
        markers(side).y = markers(side).y + head - markers(side).head;

%         if is_debug
%             figure(10+side);
%             clf;
%             hold on
%             plot3(markers(side).head(1),markers(side).head(2),markers(side).head(3)+0.1,'.','MarkerEdgeColor','r','MarkerSize',20);
%             plot3(markers(side).tail(1),markers(side).tail(2),markers(side).tail(3)+0.1,'.','MarkerEdgeColor','r','MarkerSize',20);
%             plot3(head(1), head(2), head(3)+0.1,'.','MarkerEdgeColor','g','MarkerSize',20);
%             plot3(tail(1), tail(2), tail(3)+0.1,'.','MarkerEdgeColor','g','MarkerSize',20);
%             surface('XData',xx,'YData',yy,'ZData',zz,'CData',imat,'FaceColor','texturemap','EdgeColor','none');
%             title(['Side: ', num2str(side), '; Green = Local Maxima Refine']);
%             hold off
%         end

        markers(side).head = head;
        markers(side).tail = tail;
        [~, ~, rmarkers] = ea_resolvecoords(markers, options, 1);

        if saveimg
            hf = figure(20+side);
            set(gcf,'Color',[0.1,0.1,0.1]);
            clf(hf);

            p(1) = subplot(2,2,[1,3]);
            plot3(head(1), t_head(2)-0.01, head(3),'.','MarkerEdgeColor','r','MarkerSize',20);
            hold on;
            plot3(tail(1), t_tail(2)-0.01, tail(3),'.','MarkerEdgeColor','g','MarkerSize',20);
            if is_debug % show post resolvecoords markers
                plot3(rmarkers(side).head(1),rmarkers(side).head(2)-0.01,rmarkers(side).head(3),'^','MarkerEdgeColor','r','MarkerSize',20);
                plot3(rmarkers(side).tail(1),rmarkers(side).tail(2)-0.01,rmarkers(side).tail(3),'^','MarkerEdgeColor','g','MarkerSize',20);
            end
            surface('XData',xx,'YData',yy,'ZData',zz,'CData',imat,'FaceColor','texturemap','EdgeColor','none');
            colormap(gray);
            hold off;
            axis tight;
            axis off;
            title('Lead-DBS Automated Reconstruction', 'Color', 'w');
            view(0,0);
            text(min(xx(:)), min(yy(:)), min(zz(:))*1.05, [options.subj.subjId, '; Side: ', num2str(side)], 'Color', 'w', 'FontSize', 12);

            p(2) = subplot(2,2,2);
            imagesc(tail_slice);
            hold on;
            plot(tx, ty, '.','MarkerEdgeColor','g','MarkerSize',20);
            hold off;
            axis tight;
            axis off;
            colormap(gray);

            p(3) = subplot(2,2,4);
            imagesc(head_slice);
            hold on;
            plot(hx, hy, '.','MarkerEdgeColor','r','MarkerSize',20);
            hold off;
            axis tight;
            axis off;
            colormap(gray);

            set(hf,'PaperUnits','inches','PaperPosition',[0 0 4 4], 'InvertHardCopy', 'off');
            print(hf, setBIDSEntity(options.subj.recon.recon, 'desc', ['electrode', num2str(side)], 'ext', '.jpg'), '-djpeg75', '-r300');
            if ~is_debug
                close(hf);
            end
        end
    end % for loop side iteration

    ea_save_reconstruction(coords_mm, trajectory, rmarkers, elmodel, 0, options);
    if can_export
        ea_exportfiducials(options, setBIDSEntity(options.subj.recon.recon, 'desc', 'electrodefiducials', 'ext', options.prefs.reco.exportfiducials));
    end

end

function [head, tail] = convert2mm(pidx, pidy, deltax, deltay, deltaz, xx, yy, zz, b)
    head = ones(1,3);
    tail = head;
    head(1) = pidx(1) * deltax + xx(2,1) - (xx(2,1) - xx(1,1))*(1-pidy(1)/size(b,1));
    head(2) = pidy(1) * deltay + yy(1);
    head(3) = pidy(1) * deltaz + zz(1);
    tail(1) = pidx(end) * deltax + xx(2,1) - (xx(2,1) - xx(1,1))*(1-pidy(end)/size(b,1));
    tail(2) = pidy(end) * deltay + yy(1);
    tail(3) = pidy(end) * deltaz + zz(1);
end

function [ax, ay, crop_slice, cx, cy] = refine_axial(vol, fiducial, is_debug)

    scales = spm_imatrix(vol.mat);
    xscale = scales(7)/2; % divide by 2 due to interpolation factor in ea_sample_slice()
    yscale = scales(8)/2;

    width = ceil(6/xscale); % mm width
    f_vox = vol.mat\[fiducial, 1]';
    f_vox = f_vox(1:3,:)';
    [slice,~,~,~]=ea_sample_slice(vol,'tra',width,'vox',f_vox,1);
    b = slice;
    b(b<(max(slice(:))*0.2)) = 0;


    [ys, xs] = find(b);
    mx = median(xs);
    my = median(ys);
    ax = (mx - 2*width);% * xscale;
    ay = (my - 2*width);% * yscale;

    [h, w] = size(slice);
    cropx = ceil(4/xscale);
    cropy = ceil(4/yscale);
    L = round(mx - cropx);
    R = round(mx + cropx);
    U = round(my + cropy);
    D = round(my - cropy);
    if (L < 1)
        L = 1;
    end
    if (R > w)
        R = w;
    end
    if (U > h)
        U = h;
    end
    if (D < 1)
        D = 1;
    end
    crop_slice = slice(D:U, L:R);
    cx = mx - L + 1;
    cy = my - D + 1;

    if is_debug
        f = figure(40);
        clf(f);
        imagesc(b);
        hold on;
        plot(xs, ys, '.','MarkerEdgeColor','r','MarkerSize',5);
        plot(mx, my, '.','MarkerEdgeColor','g','MarkerSize',20);
        plot(f_vox(1), f_vox(2), '.','MarkerEdgeColor','k','MarkerSize',20);
        hold off;
        f = figure(41);
        clf(f);
        imagesc(crop_slice);
        hold on;
        plot(cx, cy, '.','MarkerEdgeColor','r','MarkerSize',20);
        hold off;
    end
end

function [cx, cy] = find_centre(slice, pidx, pidy, elheight, is_debug) % elheight in voxels
    [h, w] = size(slice);
    cx = pidx;
    cy = cx;
    for i = 1:length(pidx)
        x = pidx(i);
        y = pidy(i);
        % define bounding box
        L = round(x - elheight*2);
        R = round(x + elheight*2);
        U = round(y + elheight*10);
        D = round(y - elheight*10);
        if (L < 1)
            L = 1;
        end
        if (R > w)
            R = w;
        end
        if (U > h)
            U = h;
        end
        if (D < 1)
            D = 1;
        end
        crop = slice(D:U, L:R);
        [ys, xs] = find(crop);
        mx = median(xs);
        my = median(ys);
        if is_debug
            figure(30+i);
            clf;
            hold on;
            imagesc(crop);
            colormap(gray);
            plot(xs, ys, '.','MarkerEdgeColor','r','MarkerSize',5);
            plot(mx, my, '.','MarkerEdgeColor','g','MarkerSize',30);
            hold off;
        end
        cx(i) = mx - (R-L+1)/2;
        cy(i) = my - (U-D+1)/2;
    end
end

% Function copied from ea_mancor_updatescene.m
function hdtrajectory=genhd_inside(trajectory)
    resolution=20;
    hdtrajectory(:,1)=interp1q([1:length(trajectory)]',trajectory(:,1),[1:1/resolution:length(trajectory)]');
    hdtrajectory(:,2)=interp1q([1:length(trajectory)]',trajectory(:,2),[1:1/resolution:length(trajectory)]');
    hdtrajectory(:,3)=interp1q([1:length(trajectory)]',trajectory(:,3),[1:1/resolution:length(trajectory)]');
end
