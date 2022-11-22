function ea_genflipspace
fliplrDir = [ea_space, 'fliplr', filesep];
if isempty(ea_regexpdir(fliplrDir, 'Composite\.nii\.gz$', 0))
    
    answ=questdlg('The transform to map from left to right hemisphere in MNI152NLin2009bAsym space needs to be installed or built. For reproducable results please install the transform', ...
    'Install Flip Transform', 'Install', 'Build', 'Cancel', 'Install');
    switch answ
        case 'Cancel'
            ea_error('User pressed cancel.');
        case 'Install'
            ea_checkinstall('fliplr');
            return
        otherwise
            ea_dispt('Generating nonlinear warp between hemispheres. This may take a while...');

            spacedef = ea_getspacedef;
            numModality = length(spacedef.templates);
            moving = cell(numModality, 1);
            template = cell(numModality, 1);

            ea_mkdir([ea_space, 'fliplr']);
            for t=1:numModality
                moving{t} = [fliplrDir, spacedef.templates{t}, '.nii'];
                template{t} = [ea_space, spacedef.templates{t}, '.nii'];
                copyfile(template{t}, moving{t});
                ea_flip_lr(moving{t}, moving{t});
            end

            output = [fliplrDir, 'output.nii'];

            weights = ones(numModality, 1) * 1.25;

            options.prefs = ea_prefs;

            ea_ants_nonlinear(template, moving, output, weights, '', options);

            move([fliplrDir, 'outputComposite.nii.gz'], [fliplrDir, 'Composite.nii.gz']);
            move([fliplrDir, 'outputInverseComposite.nii.gz'], [fliplrDir, 'InverseComposite.nii.gz']);

            json.method = 'ANTs (Avants 2008)';
            json.custom = 1;
            savejson('', json, [fliplrDir, 'normmethod.json']);

            ea_delete(moving);
            ea_delete(output);
            delete([fliplrDir, 'ea_ants_command.txt']);

            ea_dispt('');
    end
end
