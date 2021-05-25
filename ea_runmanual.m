function [coords_mm,trajectory,markers]=ea_runmanual(options)

directory = [options.root,options.patientname,filesep];
load([ea_getearoot,'templates',filesep,'electrode_models',filesep,options.elspec.matfname,'.mat']);

normdist=pdist([electrode.head_position;electrode.tail_position]);


for side=options.sides
    spm('defaults', 'fmri');
    Fgraph = spm_figure('GetWin', 'Graphics');
    Finter = spm('FnUIsetup','Select Trajectory', 0);

    figure(Fgraph); clf;

    switch options.modality
        case 1 % MRI
            vol = spm_vol([directory,options.prefs.tranii_unnormalized]);
        case 2 % CT
            vol = spm_vol([directory,options.prefs.ctnii_coregistered]);
    end
    spm_orthviews('Reset');
    h = spm_orthviews('Image', vol);
    colormap('gray');
    cameratoolbar('resetcamera')
    cameratoolbar('close')
    rotate3d off;

    if spm_input('Adjust contrast', 1.5, 'y/n', [1,0],2)
        pc = spm_input('Percentiles', 1.5, 'w', '3 97', 2, 100);
        wn = spm_summarise(vol, 'all', @(X) spm_percentile(X, pc));
        spm_orthviews('window', h, wn);
    end

    if spm_input('Select tip and click' , 1.5, 'OK|Retry', [1,0], 1)
        markers(side).head = spm_orthviews('Pos')';
        %spm_orthviews('Reset');
    end

    if spm_input('Select point on trajectory and click', 1.5, 'OK|Retry', [1,0], 1)
        markers(side).tail = spm_orthviews('Pos')';
        %spm_orthviews('Reset');
    end

    markers(side).tail=markers(side).head+...
        (markers(side).tail-markers(side).head)*...
        normdist/norm(...
        (markers(side).tail-markers(side).head));

    % add x and y
    [xunitv, yunitv] = ea_calcxy(markers(side).head, markers(side).tail);
    markers(side).x = markers(side).head + xunitv*(options.elspec.lead_diameter/2);
    markers(side).y = markers(side).head + yunitv*(options.elspec.lead_diameter/2);
end
[coords_mm,trajectory,markers]=ea_resolvecoords(markers,options,0);

close(Fgraph);
close(Finter);
