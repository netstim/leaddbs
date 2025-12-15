function [redo, options]=ea_prepare_dti(options)
% general DTI preprocessing goes here

redo = 0;  % Initialize return value
directory=[options.root,options.patientname,filesep];

% BIDS: If DWI files don't exist, copy from rawdata with BIDS-compliant names
if ~exist([directory, options.prefs.dti], 'file')
    options = ea_prepare_dti_bids(options);
end

if ~exist([options.root,options.patientname,filesep,options.prefs.b0],'file')
    disp('Building DTI files...');
    try %unring
        dti=ea_load_untouch_nii([options.root,options.patientname,filesep,options.prefs.dti]);
        dti.img=ea_unring(dti.img);
        ea_save_untouch_nii(dti,[options.root,options.patientname,filesep,options.prefs.dti]);
    end

    % export B0
    if ~exist([options.root,options.patientname,filesep,options.prefs.b0],'file')
        ea_exportb0(options);
    else
        % BIDS FIX: Check and correct b0 header even if file exists
        % This fixes issues with DWI data that have incorrect qform/sform origins
        b0Path = [options.root,options.patientname,filesep,options.prefs.b0];
        try
            V = spm_vol(b0Path);
            origin_vox = [-V.mat(1,4)/V.mat(1,1), -V.mat(2,4)/V.mat(6), -V.mat(3,4)/V.mat(11)];
            center_vox = V.dim / 2;
            
            % If origin is more than 20% off from center, correct it
            if abs(origin_vox(3) - center_vox(3)) > 0.2 * center_vox(3)
                fprintf('\nCorrecting existing b0 header origin (was at voxel [%.1f, %.1f, %.1f], center is [%.1f, %.1f, %.1f])...\n', ...
                    origin_vox(1), origin_vox(2), origin_vox(3), center_vox(1), center_vox(2), center_vox(3));
                
                % Read the image
                img_data = spm_read_vols(V);
                
                % Create new header with corrected origin (centered)
                V_new = V;
                V_new.mat(1,4) = -V.mat(1,1) * center_vox(1);
                V_new.mat(2,4) = -V.mat(6) * center_vox(2);
                V_new.mat(3,4) = -V.mat(11) * center_vox(3);
                
                % Write with corrected header
                spm_write_vol(V_new, img_data);
                fprintf('Existing b0 header corrected successfully.\n\n');
            end
        catch ME
            warning('Failed to check/correct existing b0 header: %s', ME.message);
        end
    end

    % export FA
    if ~exist([options.root,options.patientname,filesep,options.prefs.fa],'file')
        ea_isolate_fa(options);
    end
else
    disp('B0 found, no need to rebuild.');
end

% restore dti
if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],'file')
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.dti)],...
        [options.root,options.patientname,filesep,options.prefs.dti]);
    redo=1; % apparently prior run crashed - redo to be safe.
end
% restore b0
if exist([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],'file')
    movefile([options.root,options.patientname,filesep,ea_stripext(options.prefs.b0)],...
        [options.root,options.patientname,filesep,options.prefs.b0]);
    redo=1; % apparently prior run crashed - redo to be safe.
end
