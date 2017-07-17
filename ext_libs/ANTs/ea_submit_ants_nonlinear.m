function ea_submit_ants_nonlinear(props)

if props.stagesep
    ea_antsnl_multistep(props) % run each ANTs stage separately
else
    ea_antsnl_monostep(props) % run all ANTs stages together
end


function ea_antsnl_monostep(props)

cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --output [',ea_path_helper(props.outputbase), ',', props.outputimage, ']', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 0', ...
    ' --float 1',...
    ' --winsorize-image-intensities [0.005,0.995]', ...
    ' --write-composite-transform 1', ...
    props.rigidstage, props.affinestage, props.synstage, props.slabstage, props.synmaskstage];

fid=fopen([props.directory,'ea_ants_command.txt'],'a');
fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end

ea_conv_antswarps(props.directory);

function ea_antsnl_multistep(props)


%% Execute stages one for one

% linear:
cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    ' --output ',ea_path_helper(props.outputbase),'', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 0', ...
    ' --winsorize-image-intensities[0.005,0.995]', ...
    ' --float 1',...
    props.rigidstage, props.affinestage];

fid=fopen([props.directory,'ea_ants_command.txt'],'a');
fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end
stack=[' --initial-moving-transform ',ea_path_helper(props.outputbase),'0GenericAffine.mat'];
tstack={[' --transform [',ea_path_helper(props.outputbase),'0GenericAffine.mat,0]']};
itstack={[' --transform [',ea_path_helper(props.outputbase),'0GenericAffine.mat,1]']};

% SyN:
cmd = [props.ANTS, ' --verbose 1', ...
    ' --dimensionality 3', ...
    stack ... % Apply Linear
    ' --output ',ea_path_helper(props.outputbase), 'Diff', ...
    ' --interpolation Linear', ...
    ' --use-histogram-matching 0', ...
    ' --winsorize-image-intensities[0.005,0.995]', ...
    ' --float 1',...
    props.synstage];

display(cmd)
fid=fopen([props.directory,'ea_ants_command.txt'],'a');
fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
fclose(fid);

if ~ispc
    system(['bash -c "', cmd, '"']);
else
    system(cmd);
end


stack=[' --initial-moving-transform ',ea_path_helper(props.outputbase),'Diff1Warp.nii.gz', ...
    ' --initial-moving-transform ',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat'];
tstack={[' --transform ',ea_path_helper(props.outputbase),'Diff1Warp.nii.gz'],...
    [' --transform [',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat,0]']};
itstack={[' --transform [',ea_path_helper(props.outputbase),'Diff0GenericAffine.mat,1]'],...
    [' --transform ',ea_path_helper(props.outputbase),'Diff1InverseWarp.nii.gz']};

% Slab:
if ~isempty(props.slabstage)
    cmd = [props.ANTS, ' --verbose 1', ...
        ' --dimensionality 3', ...
        stack ... % Apply Syn, Apply Linear
        ' --output ',ea_path_helper(props.outputbase), 'DiffSlab', ...
        ' --interpolation Linear', ...
        ' --use-histogram-matching 0', ...
        ' --winsorize-image-intensities[0.005,0.995]', ...
        ' --float 1',...
        props.slabstage];
    
    display(cmd)
    fid=fopen([props.directory,'ea_ants_command.txt'],'a');
    fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
    fclose(fid);
    
    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
    stack=[' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab2Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab1Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat'];
    tstack={[' --transform ',ea_path_helper(props.outputbase),'DiffSlab2Warp.nii.gz'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffSlab1Warp.nii.gz'],...
        [' --transform [',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat,0]']};
    itstack={[' --transform [',ea_path_helper(props.outputbase),'DiffSlab0GenericAffine.mat,1]'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffSlab2InverseWarp.nii.gz']};
end


% SyNSubcorticalRefineStage:
if ~isempty(props.synmaskstage)
    cmd = [props.ANTS, ' --verbose 1', ...
        ' --dimensionality 3', ...
        stack ... % Apply Slab % Apply Syn % Apply Linear
        ' --output ',ea_path_helper(props.outputbase), 'DiffScrf', ...
        ' --interpolation Linear', ...
        ' --use-histogram-matching 0', ...
        ' --winsorize-image-intensities[0.005,0.995]', ...
        ' --float 1',...
        props.synmaskstage];
    
    display(cmd)
    fid=fopen([props.directory,'ea_ants_command.txt'],'a');
    fprintf(fid,[datestr(datetime('now')),':\n',cmd,'\n\n']);
    fclose(fid);
    
    if ~ispc
        system(['bash -c "', cmd, '"']);
    else
        system(cmd);
    end
    stack=[' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf2Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf1Warp.nii.gz',...
        ' --initial-moving-transform ',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat'];
    tstack={[' --transform ',ea_path_helper(props.outputbase),'DiffScrf2Warp.nii.gz'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffScrf1Warp.nii.gz'],...
        [' --transform [',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat,0]']};
    itstack={[' --transform [',ea_path_helper(props.outputbase),'DiffScrf0GenericAffine.mat,1]'],...
        [' --transform ',ea_path_helper(props.outputbase),'DiffScrf2InverseWarp.nii.gz']};
end

%% combine all warps to a single one:
if ispc
    sufx='.exe';
else
    sufx=computer('arch');
end
antsApply=[ea_getearoot,'ext_libs',filesep,'ANTS',filesep,'antsApplyTransforms.',sufx];
template=[ea_space,'t1.nii'];
outputformat='.nii.gz';

tstring=[]; itstring=[];
for t=1:length(tstack)-1
    tstring=[tstring,tstack{t}];
    itstring=[itstring,itstack{t}];
end
options=ea_getptopts(props.directory);
prenii=[props.directory,options.prefs.prenii_unnormalized];

% apply command:
cmd=[antsApply,' -r ',template,...
    tstring, ...
    ' -o [',ea_path_helper([props.directory,'glanatComposite',outputformat]),',1]'];
icmd=[antsApply,' -r ',ea_path_helper(prenii),...
    tstring, ...
    ' -o [',ea_path_helper([props.directory,'glanatInverseComposite',outputformat]),',1]'];

if exist('cmd','var')
    if ~ispc
        system(['bash -c "', cmd, '"']);
        system(['bash -c "', icmd, '"']);
    else
        system(cmd);
        system(icmd);
    end
end

% delete all old-version warps
warning('off')
ea_delete([props.directory,'glanatComposite.h5']);
ea_delete([props.directory,'glanatInverseComposite.h5']);
ea_delete([props.directory,'glanat0GenericeAffine.mat']);
ea_delete([props.directory,'glanat*Warp.mat']);
ea_delete([props.directory,'glanat*InverseWarp.mat']);
warning('on')









