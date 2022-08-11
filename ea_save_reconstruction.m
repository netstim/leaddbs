function ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,manually_corrected,options,outFilePath)

if ~exist('outfname', 'var')
    outFilePath = options.subj.recon.recon;
end

ea_mkdir(fileparts(outFilePath));

for side=options.sides
    reco.props(side).elmodel=elmodel;
    reco.props(side).manually_corrected=manually_corrected;
end

if options.native
    reco.native.coords_mm=coords_mm;
    reco.native.trajectory=trajectory;
    reco.native.markers=markers;
    checkxy = 1;
    if checkxy % checks for problems with the length of markers x and y - could be removed at some point
        for numleads = 1:length(reco.native.markers)
            x =  reco.native.markers(numleads).x - reco.native.markers(numleads).head;
            y =  reco.native.markers(numleads).y - reco.native.markers(numleads).head;
            x = (x./norm(x)) .* (options.elspec.lead_diameter/2);
            y = (y./norm(y)) .* (options.elspec.lead_diameter/2);
            reco.native.markers(numleads).x = reco.native.markers(numleads).head + x;
            reco.native.markers(numleads).y = reco.native.markers(numleads).head + y;
            clear x y
        end
    end
    save(outFilePath, 'reco');

    if isfield(options,'hybridsave')
        ea_dispt('Warping fiducials to template space');
        ea_reconstruction2mni(options);
        if options.prefs.reco.saveACPC
            ea_dispt('Mapping fiducials to AC/PC space');
            ea_reconstruction2acpc(options);
            ea_dispt('');
        end
        load(outFilePath, 'reco');
        [reco,corrected]=ea_checkswap_lr(reco,options); % PaCER support, right could be left and vice versa.

        save(outFilePath, 'reco');
        ea_dispt('');
    end
else
    reco.mni.coords_mm=coords_mm;
    reco.mni.trajectory=trajectory;
    reco.mni.markers=markers;
    save(outFilePath, 'reco');

    if isfield(options,'hybridsave')
        ea_dispt('Warping fiducials to native space');
        ea_reconstruction2native(options);
        if options.prefs.reco.saveACPC
            ea_dispt('Mapping fiducials to AC/PC space');
            ea_reconstruction2acpc(options);
            ea_dispt('');
        end
        load(outFilePath, 'reco');
        [reco,corrected]=ea_checkswap_lr(reco,options); % PaCER support, right could be left and vice versa.

        save(outFilePath, 'reco');
        ea_dispt('');

        if corrected
            options.hybridsave=1;
            ea_save_reconstruction(reco.mni.coords_mm,reco.mni.trajectory,reco.mni.markers,elmodel,manually_corrected,options)
        end
    end
end


function [reco,corrected]=ea_checkswap_lr(reco,options)
options.native=0; % this can only be done in MNI space.
corrected=0;
if length(reco.mni.coords_mm)==2 && ~any(cellfun(@isempty,reco.mni.coords_mm))
    if mean(reco.mni.coords_mm{1}(:,1))<mean(reco.mni.coords_mm{2}(:,1)) % RL swapped
        % swap RL:
        options.hybridsave=1;
        ncoords_mm{1}=reco.mni.coords_mm{2};
        ncoords_mm{2}=reco.mni.coords_mm{1};
        ntrajectory{1}=reco.mni.trajectory{2};
        ntrajectory{2}=reco.mni.trajectory{1};
        nmarkers(1)=reco.mni.markers(2);
        nmarkers(2)=reco.mni.markers(1);

        reco.mni.coords_mm=ncoords_mm;
        reco.mni.trajectory=ntrajectory;
        reco.mni.markers=nmarkers;

        ncoords_mm{1}=reco.native.coords_mm{2};
        ncoords_mm{2}=reco.native.coords_mm{1};
        ntrajectory{1}=reco.native.trajectory{2};
        ntrajectory{2}=reco.native.trajectory{1};
        nmarkers(1)=reco.native.markers(2);
        nmarkers(2)=reco.native.markers(1);

        reco.native.coords_mm=ncoords_mm;
        reco.native.trajectory=ntrajectory;
        reco.native.markers=nmarkers;

        corrected=1;
    end
end

vizz=0;

% check that markers are correct (important for directional leads):
if ~reco.props(options.sides(1)).manually_corrected
    options.hybridsave=1;

    for side=options.sides
        if vizz
            figure
            hold on
            plot3(reco.mni.trajectory{side}(:,1),reco.mni.trajectory{side}(:,2),reco.mni.trajectory{side}(:,3),'r-');
            plot3(reco.mni.markers(side).head(:,1),reco.mni.markers(side).head(:,2),reco.mni.markers(side).head(:,3),'y*');
            plot3(reco.mni.markers(side).tail(:,1),reco.mni.markers(side).tail(:,2),reco.mni.markers(side).tail(:,3),'m*');
            plot3(reco.mni.markers(side).x(:,1),reco.mni.markers(side).x(:,2),reco.mni.markers(side).x(:,3),'k*');
            plot3(reco.mni.markers(side).y(:,1),reco.mni.markers(side).y(:,2),reco.mni.markers(side).y(:,3),'g*');
        end

        if reco.mni.markers(side).head(2)>reco.mni.markers(side).y(2) % Flip y marker if head is anterior to y marker.
            reco.mni.markers(side).y=reco.mni.markers(side).head+(reco.mni.markers(side).head-reco.mni.markers(side).y); % 180 deg flip
            corrected=1;
        end

        if reco.mni.markers(side).head(1)>reco.mni.markers(side).x(1) % Flip x marker if head is right to x marker.
            reco.mni.markers(side).x=reco.mni.markers(side).head+(reco.mni.markers(side).head-reco.mni.markers(side).x);
            corrected=1;
        end

        if vizz
            plot3(reco.mni.markers(side).x(:,1),reco.mni.markers(side).x(:,2),reco.mni.markers(side).x(:,3),'ko');
            plot3(reco.mni.markers(side).y(:,1),reco.mni.markers(side).y(:,2),reco.mni.markers(side).y(:,3),'go');
            axis equal
            keyboard
        end
    end
    reco.mni.coords_mm=ea_resolvecoords(reco.mni.markers,options);
end
