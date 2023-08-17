function varargout =  ea_runslicer(options, task)
%% Function to launch slicer and load *.nii files
%  Last Revision: 21/05/2018
%  Thushara Perera (c) 2018 Bionics Institute
%  Input:
%   - lead dbs options struct
%   - task integer ID to tell function which files to open
%     (see inline comments for more detail)
%  Output:
%   - relevant slicer scene (mrml) file will be saved in patient folder
%   - slicer fiducial marker file (fcsv) will be saved in patient folder
%   - temporary python script file will be written to lead_dbs root
%
%  Things to note: user must set the path of the slicer executable. It would
%  be nice if there was an easy method to determine the path automatically.

    options.prefs = ea_prefs('');
    if ~isfield(options.prefs,'slicer') || isempty(options.prefs.slicer.dir)
        ea_hastoolbox('slicer');
        options.prefs = ea_prefs('');
        if ~isfield(options.prefs, 'slicer')
            warning(sprintf('3D Slicer path not set!\nPlease set ''prefs.slicer.dir'' in your preference file.'))
            return;
        end
    end

    slicer_path = options.prefs.slicer.dir;
    if isempty(slicer_path)
        warning(sprintf('3D Slicer path not set!\nPlease set ''prefs.slicer.dir'' in your preference file.'))
        return;
    end

    % 'prefs.slicer.dir' can be either the folder containing Slicer executable or
    % the full path to the excutable itself
    if ismac
        if isfolder(slicer_path) && regexp(slicer_path, 'Slicer\.app/?$')
            SLICER = fullfile(slicer_path,'Contents','MacOS','Slicer');
        elseif isfile(slicer_path) && regexp(slicer_path, 'Slicer$')
            SLICER = slicer_path;
        else
            warning('Slicer executable not found!');
            return;
        end
    elseif isunix
        if isfolder(slicer_path) && isfile(fullfile(slicer_path, 'Slicer'))
            SLICER = fullfile(slicer_path, 'Slicer');
        elseif isfile(slicer_path) && regexp(slicer_path, 'Slicer$')
            SLICER = slicer_path;
        else
            warning('Slicer executable not found!');
            return;
        end
    elseif ispc
        if isfolder(slicer_path) && isfile(fullfile(slicer_path, 'Slicer.exe'))
            SLICER = fullfile(slicer_path, 'Slicer.exe');
        elseif isfile(slicer_path) && regexp(slicer_path, 'Slicer\.exe$')
            SLICER = slicer_path;
        else
            warning('Slicer executable not found!');
            return;
        end
    end

    if isempty(options.uipatdirs)
        warning('No patient selected!');
        return;
    elseif length(options.uipatdirs)>1 && task ~= 5
        warning('Slicer module only works for single patient!');
        return;
    end

    switch task
        case -1 % launch slicer for lead reconstruction
            mrmltag = 'reconstruction';
            allfiles = struct2cell(options.subj.coreg.anat.postop);
            allfiles = allfiles(~contains(allfiles, 'tonemapped'));

            filepaths = allfiles(isfile(allfiles));
            filenames = regexp(filepaths, ['(?<=\', filesep, ')[\w-]+(?=\.nii(\.gz)?$)'], 'match', 'once');
            nfiles = length(filenames);

            WriteReconstructionFiducialFile(options);

        case 1 % show original volumes
            mrmltag = 'original';
            allfiles = [struct2cell(options.subj.preproc.anat.postop)
                struct2cell(options.subj.preproc.anat.preop)];

            filepaths = allfiles(isfile(allfiles));
            filenames = regexp(filepaths, ['(?<=\', filesep, ')[\w-]+(?=\.nii(\.gz)?$)'], 'match', 'once');
            nfiles = length(filenames);

        case 2 % show data after co-registration
            mrmltag = 'coregistered';
            allfiles = [struct2cell(options.subj.coreg.anat.postop)
                struct2cell(options.subj.coreg.anat.preop)];
            allfiles = allfiles(~contains(allfiles, 'tonemapped'));

            filepaths = allfiles(isfile(allfiles));
            filenames = regexp(filepaths, ['(?<=\', filesep, ')[\w-]+(?=\.nii(\.gz)?$)'], 'match', 'once');
            nfiles = length(filenames);

        case 3 % show data after normalization
            mrmltag = 'normalized';
            allfiles = [[ea_space, 't2.nii']
                struct2cell(options.subj.norm.anat.postop)
                struct2cell(options.subj.norm.anat.preop)];
            allfiles = allfiles(~contains(allfiles, 'tonemapped'));

            filepaths = allfiles(isfile(allfiles));
            filenames = regexp(filepaths, ['(?<=\', filesep, ')[\w-]+(?=\.nii(\.gz)?$)'], 'match', 'once');
            nfiles = length(filenames);

        case 4 % show electrode localization
            if isfile(options.subj.recon.recon)
                options.native = 0; % Export fiducial only in MNI space
                ea_exportfiducials(options, setBIDSEntity(options.subj.recon.recon, 'desc', 'electrodefiducials', 'ext', 'fcsv'));
            else
                warning('Please run reconstruction first...');
                return;
            end
            mrmltag = 'electrodes';
            allfiles = [[ea_space, 't2.nii']
                struct2cell(options.subj.norm.anat.postop)
                struct2cell(options.subj.norm.anat.preop)];
            allfiles = allfiles(~contains(allfiles, 'tonemapped'));

            filepaths = allfiles(isfile(allfiles));
            filenames = regexp(filepaths, ['(?<=\', filesep, ')[\w-]+(?=\.nii(\.gz)?$)'], 'match', 'once');
            nfiles = length(filenames);

        case 5 % return slicer path
            varargout = {SLICER};
            return

        otherwise
            warning('Task ID not recognised');
            return;
    end

    slicer_folder = [options.subj.subjDir, filesep, 'slicer'];
    ea_mkdir(slicer_folder);
    script_path = [slicer_folder, filesep, 'slicer.py'];
    scene_path = [slicer_folder, filesep, 'sub-', options.subj.subjId, '_desc-', mrmltag, '.mrml'];
    if (nfiles < 1)
        warning('Need at least one volume image to load into slicer');
        return;
    end

    fid = fopen(scene_path, 'w');
    if (task == 4 || task == -1)
        fprintf(fid, [GetFiducialBeginning(), '\r\n\r\n']);
    else
        fprintf(fid, [GetBeginning(), '\r\n\r\n']);
    end
    for i=1:nfiles
        fprintf(fid, [GetFileXML(i, filepaths{i}, filenames{i}), '\r\n\r\n']);
    end
    if (task == 4 || task == -1)
        fprintf(fid, GetFiducialEnding(setBIDSEntity(options.subj.recon.recon, 'desc', 'electrodefiducials', 'ext', 'fcsv')));
    else
        fprintf(fid, GetEnding());
    end
    fclose(fid);

    fid = fopen(script_path, 'w'); % write temp python script to load volumes
    fprintf(fid, ['slicer.util.loadScene("', strrep(scene_path, '\', '/'), '")\r\n']);
    fclose(fid);
    disp('Loading up 3D Slicer...');
    ea_libs_helper('', 'unset');
    if (task > 0)
        system(['"', SLICER, '" --no-splash --python-script "', script_path, '" &']);
        % the trailing '&' returns control back to matlab without waiting for slicer to close
    else
        system(['"', SLICER, '" --no-splash --python-script "', script_path, '"']);
    end

    % ea_delete(script_path);
    % Please do not delete the above script. Slicer runs asynchronously in
    % the background. A race condition will develop and the script will be
    % deleted by Matlab before Slicer has a chance to read it!
end

function WriteReconstructionFiducialFile(options)
    fiducial_path = setBIDSEntity(options.subj.recon.recon, 'desc', 'electrodefiducials', 'ext', 'fcsv');
    header = ['# Markups fiducial file version = 4.7\r\n',...
              '# CoordinateSystem = 0\r\n',...
              '# columns = id,x,y,z,ow,ox,oy,oz,vis,sel,lock,label,desc,associatedNodeID\r\n'];
    fid = fopen(fiducial_path, 'w');
    fprintf(fid, header);

    counter = 0;
    if isfile(options.subj.recon.recon)
        options.native = 1;
        options.loadnativereco = 1; % Load native reco intead of scrf
        [~,~,markers] = ea_load_reconstruction(options);
        for side = options.sides
            h = markers(side).head;
            t = markers(side).tail;
            fprintf(fid, ['vtkMRMLMarkupsFiducialNode_',num2str(counter),',',num2str(h(1)),',',num2str(h(2)),',',num2str(h(3)),...
                ',0,0,0,1,1,1,0,Head_',num2str(side),',,vtkMRMLScalarVolumeNode2\r\n']);
            fprintf(fid, ['vtkMRMLMarkupsFiducialNode_',num2str(counter+1),',',num2str(t(1)),',',num2str(t(2)),',',num2str(t(3)),...
                ',0,0,0,1,1,1,0,Tail_',num2str(side),',,vtkMRMLScalarVolumeNode2\r\n']);
            counter = counter + 2;
        end
    else
        h = 14; t = 9;
        for side = options.sides
            fprintf(fid, ['vtkMRMLMarkupsFiducialNode_',num2str(counter),',',num2str(h),',-11,-5',...
                ',0,0,0,1,1,1,0,Head_',num2str(side),',,vtkMRMLScalarVolumeNode2\r\n']);
            fprintf(fid, ['vtkMRMLMarkupsFiducialNode_',num2str(counter+1),',',num2str(t),',-15,-15',...
                ',0,0,0,1,1,1,0,Tail_',num2str(side),',,vtkMRMLScalarVolumeNode2\r\n']);
            counter = counter + 2;
            if mod(side,2) == 0 % is even
                h = h + 10;
                t = t + 10;
            end
            % flip left/right after each iteration
            h = -h;
            t = -t;
        end
    end

    fclose(fid);
end


function txt = GetFileXML(index, filepath, name)

    path = strrep(filepath, '\', '\\'); % avoid escape char errors

    vdisplay = ['<VolumeDisplay id="vtkMRMLScalarVolumeDisplayNode', num2str(index), '" ',...
        'name="VolumeDisplay" hideFromEditors="true" selectable="true" selected="false" color="0.5 0.5 0.5" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" edgeVisibility="false" clipping="false" sliceIntersectionVisibility="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="true" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" colorNodeID="vtkMRMLColorTableNodeGrey"  window="100" level="50" upperThreshold="32767" lowerThreshold="-32768" interpolate="1" autoWindowLevel="0" applyThreshold="0" autoThreshold="0" ></VolumeDisplay>'];


    vstorage = ['<VolumeArchetypeStorage ',...
        'id="vtkMRMLVolumeArchetypeStorageNode', num2str(index),'" ',...
        'name="VolumeArchetypeStorage_', num2str(index),'" ',...
        'hideFromEditors="true" selectable="true" selected="false" ',...
        'fileName="',path,'" ',...
        'useCompression="1" defaultWriteFileExtension="nrrd" readState="0" writeState="0" centerImage="0" UseOrientationFromFile="1" ></VolumeArchetypeStorage>'];


    volume = ['<Volume id="vtkMRMLScalarVolumeNode', num2str(index),'" ',...
        'name="',name,'" ',...
        'hideFromEditors="false" selectable="true" selected="false" ',...
        'displayNodeRef="vtkMRMLScalarVolumeDisplayNode', num2str(index),'" ',...
        'storageNodeRef="vtkMRMLVolumeArchetypeStorageNode', num2str(index),'" ',...
        'references="display:vtkMRMLScalarVolumeDisplayNode', num2str(index),';',...
        'storage:vtkMRMLVolumeArchetypeStorageNode', num2str(index),';" ',...
        'userTags="" ></Volume>'];

    txt = [vdisplay, vstorage, volume];

end


function txt = GetBeginning()
     txt = ['<MRML  version="Slicer4.4.0" userTags="">',...
     '<Crosshair id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="JumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"></Crosshair>',...
     '<Selection id="vtkMRMLSelectionNodeSingleton" name="Selection" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" frequencyUnitNodeRef="vtkMRMLUnitNodeApplicationFrequency" intensityUnitNodeRef="vtkMRMLUnitNodeApplicationIntensity" lengthUnitNodeRef="vtkMRMLUnitNodeApplicationLength" timeUnitNodeRef="vtkMRMLUnitNodeApplicationTime" velocityUnitNodeRef="vtkMRMLUnitNodeApplicationVelocity" references="unit/frequency:vtkMRMLUnitNodeApplicationFrequency;unit/intensity:vtkMRMLUnitNodeApplicationIntensity;unit/length:vtkMRMLUnitNodeApplicationLength;unit/time:vtkMRMLUnitNodeApplicationTime;unit/velocity:vtkMRMLUnitNodeApplicationVelocity;" activeVolumeID="vtkMRMLScalarVolumeNode2" secondaryVolumeID="NULL" activeLabelVolumeID="NULL" activeFiducialListID="NULL" activePlaceNodeID="NULL" activePlaceNodeClassName="NULL" activeROIListID="NULL" activeCameraID="NULL" activeTableID="NULL" activeViewID="NULL" activeLayoutID="NULL" ></Selection>',...
     '<Interaction id="vtkMRMLInteractionNodeSingleton" name="Interaction" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" currentInteractionMode="ViewTransform" placeModePersistence="false" lastInteractionMode="ViewTransform" ></Interaction>',...
     '<View id="vtkMRMLViewNode1" name="View1" hideFromEditors="false" selectable="true" selected="false" singletonTag="1" attributes="MappedInLayout:1" layoutLabel="1" layoutName="1" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="200" letterSize="0.05" boxVisible="true" fiducialsVisible="true" fiducialLabelsVisible="true" axisLabelsVisible="true" axisLabelsCameraDependent="true" animationMode="Off" viewAxisMode="LookFrom" spinDegrees="2" spinMs="5" spinDirection="YawLeft" rotateDegrees="5" rockLength="200" rockCount="0" stereoType="NoStereo" renderMode="Perspective" useDepthPeeling="0" ></View>',...
     '<Slice id="vtkMRMLSliceNodeRed" name="Red" hideFromEditors="false" selectable="true" selected="false" singletonTag="Red" attributes="MappedInLayout:1" layoutLabel="R" layoutName="Red" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="318.002 259.289 0.5" dimensions="1002 817 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="318.002 259.289 0.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 1 0 -1.77324e-006 0 0 1 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.952941 0.290196 0.2" orientation="Axial" orientationReference="Axial" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>',...
     '<Slice id="vtkMRMLSliceNodeYellow" name="Yellow" hideFromEditors="false" selectable="true" selected="false" singletonTag="Yellow" attributes="MappedInLayout:1" layoutLabel="Y" layoutName="Yellow" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="0 0 1 0.218502 -1 0 0 -1.77324e-006 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.929412 0.835294 0.298039" orientation="Sagittal" orientationReference="Sagittal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>',...
     '<Slice id="vtkMRMLSliceNodeGreen" name="Green" hideFromEditors="false" selectable="true" selected="false" singletonTag="Green" attributes="MappedInLayout:1" layoutLabel="G" layoutName="Green" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 0 1 0.218498 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.431373 0.690196 0.294118" orientation="Coronal" orientationReference="Coronal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>',...
     '<Layout id="vtkMRMLLayoutNodevtkMRMLLayoutNode" name="Layout" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLLayoutNode" currentViewArrangement="3" guiPanelVisibility="1" bottomPanelVisibility ="1" guiPanelLR="0" collapseSliceControllers="0"',...
     ' numberOfCompareViewRows="1" numberOfCompareViewColumns="1" numberOfLightboxRows="6" numberOfLightboxColumns="6" mainPanelSize="400" secondaryPanelSize="400" ></Layout>',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeRed" name="SliceComposite" hideFromEditors="true" selectable="true" selected="false" singletonTag="Red" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Red" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeYellow" name="SliceComposite_1" hideFromEditors="true" selectable="true" selected="false" singletonTag="Yellow" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Yellow" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeGreen" name="SliceComposite_2" hideFromEditors="true" selectable="true" selected="false" singletonTag="Green" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Green" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>',...
     '<Camera id="vtkMRMLCameraNode1" name="Default Scene Camera" hideFromEditors="false" selectable="true" selected="false" userTags="" position="-174.789 463.886 65.26" focalPoint="0 0 0" viewUp="0.0301267 -0.128106 0.991303" parallelProjection="false" parallelScale="1" activetag="vtkMRMLViewNode1" appliedTransform="1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1" ></Camera>',...
     '<SubjectHierarchy id="vtkMRMLSubjectHierarchyNode1" name="SubjectHierarchy" hideFromEditors="false" selectable="true" selected="false" attributes="SubjectHierarchyVersion:2" >',...
       '<SubjectHierarchyItem id="3" name="Scene" parent="0" type="" expanded="true" attributes="Level^Scene|">',...
         '<SubjectHierarchyItem id="9" dataNode="vtkMRMLScalarVolumeNode1" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>',...
         '<SubjectHierarchyItem id="12" dataNode="vtkMRMLScalarVolumeNode2" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>',...
         '<SubjectHierarchyItem id="13" dataNode="vtkMRMLSceneViewNode1" parent="3" type="SceneViews" expanded="true"></SubjectHierarchyItem></SubjectHierarchyItem></SubjectHierarchy>',...
     '<ClipModels id="vtkMRMLClipModelsNodevtkMRMLClipModelsNode" name="ClipModels" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLClipModelsNode" clipType="0" redSliceClipState="0" yellowSliceClipState="0" greenSliceClipState="0" ></ClipModels>',...
     '<ScriptedModule id="vtkMRMLScriptedModuleNodeDataProbe" name="ScriptedModule" hideFromEditors="true" selectable="true" selected="false" singletonTag="DataProbe" ModuleName ="DataProbe" ></ScriptedModule>'...
     ];
end


function txt = GetEnding()
    txt = ['<Crosshair id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="JumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"  ></Crosshair>',...
    '<Selection id="vtkMRMLSelectionNodeSingleton" name="Selection" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" frequencyUnitNodeRef="vtkMRMLUnitNodeApplicationFrequency" intensityUnitNodeRef="vtkMRMLUnitNodeApplicationIntensity" lengthUnitNodeRef="vtkMRMLUnitNodeApplicationLength" timeUnitNodeRef="vtkMRMLUnitNodeApplicationTime" velocityUnitNodeRef="vtkMRMLUnitNodeApplicationVelocity" references="unit/frequency:vtkMRMLUnitNodeApplicationFrequency;unit/intensity:vtkMRMLUnitNodeApplicationIntensity;unit/length:vtkMRMLUnitNodeApplicationLength;unit/time:vtkMRMLUnitNodeApplicationTime;unit/velocity:vtkMRMLUnitNodeApplicationVelocity;" activeVolumeID="vtkMRMLScalarVolumeNode2" secondaryVolumeID="vtkMRMLScalarVolumeNode2" activeLabelVolumeID="NULL" activeFiducialListID="NULL" activePlaceNodeID="NULL" activePlaceNodeClassName="NULL" activeROIListID="NULL" activeCameraID="NULL" activeTableID="NULL" activeViewID="NULL" activeLayoutID="NULL"  ></Selection>',...
    '<Interaction id="vtkMRMLInteractionNodeSingleton" name="Interaction" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" currentInteractionMode="ViewTransform" placeModePersistence="false" lastInteractionMode="ViewTransform"  ></Interaction>',...
    '<View id="vtkMRMLViewNode1" name="View1" hideFromEditors="false" selectable="true" selected="false" singletonTag="1" attributes="MappedInLayout:1" layoutLabel="1" layoutName="1" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="200" letterSize="0.05" boxVisible="true" fiducialsVisible="true" fiducialLabelsVisible="true" axisLabelsVisible="true" axisLabelsCameraDependent="true" animationMode="Off" viewAxisMode="LookFrom" spinDegrees="2" spinMs="5" spinDirection="YawLeft" rotateDegrees="5" rockLength="200" rockCount="0" stereoType="NoStereo" renderMode="Perspective" useDepthPeeling="0"  ></View>',...
    '<Slice id="vtkMRMLSliceNodeRed" name="Red" hideFromEditors="false" selectable="true" selected="false" singletonTag="Red" attributes="MappedInLayout:1" layoutLabel="R" layoutName="Red" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="318.002 259.289 0.5" dimensions="1002 817 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="318.002 259.289 0.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 1 0 -1.77324e-006 0 0 1 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.952941 0.290196 0.2" orientation="Axial" orientationReference="Axial" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Slice id="vtkMRMLSliceNodeYellow" name="Yellow" hideFromEditors="false" selectable="true" selected="false" singletonTag="Yellow" attributes="MappedInLayout:1" layoutLabel="Y" layoutName="Yellow" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="0 0 1 0.218502 -1 0 0 -1.77324e-006 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.929412 0.835294 0.298039" orientation="Sagittal" orientationReference="Sagittal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Slice id="vtkMRMLSliceNodeGreen" name="Green" hideFromEditors="false" selectable="true" selected="false" singletonTag="Green" attributes="MappedInLayout:1" layoutLabel="G" layoutName="Green" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 0 1 0.218498 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.431373 0.690196 0.294118" orientation="Coronal" orientationReference="Coronal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Layout id="vtkMRMLLayoutNodevtkMRMLLayoutNode" name="Layout" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLLayoutNode" currentViewArrangement="3" guiPanelVisibility="1" bottomPanelVisibility ="1" guiPanelLR="0" collapseSliceControllers="0"',...
    ' numberOfCompareViewRows="1" numberOfCompareViewColumns="1" numberOfLightboxRows="6" numberOfLightboxColumns="6" mainPanelSize="400" secondaryPanelSize="400"  ></Layout>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeRed" name="SliceComposite" hideFromEditors="true" selectable="true" selected="false" singletonTag="Red" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Red" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeYellow" name="SliceComposite_1" hideFromEditors="true" selectable="true" selected="false" singletonTag="Yellow" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Yellow" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeGreen" name="SliceComposite_2" hideFromEditors="true" selectable="true" selected="false" singletonTag="Green" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Green" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<Camera id="vtkMRMLCameraNode1" name="Default Scene Camera" hideFromEditors="false" selectable="true" selected="false" userTags="" position="-174.789 463.886 65.26" focalPoint="0 0 0" viewUp="0.0301267 -0.128106 0.991303" parallelProjection="false" parallelScale="1" activetag="vtkMRMLViewNode1" appliedTransform="1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1"  ></Camera>',...
    '<ClipModels id="vtkMRMLClipModelsNodevtkMRMLClipModelsNode" name="ClipModels" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLClipModelsNode" clipType="0" redSliceClipState="0" yellowSliceClipState="0" greenSliceClipState="0"  ></ClipModels>',...
    '<ScriptedModule id="vtkMRMLScriptedModuleNodeDataProbe" name="ScriptedModule" hideFromEditors="true" selectable="true" selected="false" singletonTag="DataProbe" ModuleName ="DataProbe"  ></ScriptedModule>',...
    '<SceneViewStorage id="vtkMRMLSceneViewStorageNode1" name="SceneViewStorage" hideFromEditors="true" selectable="true" selected="false" fileName="Master Scene View.png" useCompression="1" defaultWriteFileExtension="png" readState="0" writeState="4" ></SceneViewStorage>',...
    '</MRML>'
    ];
end

function txt = GetFiducialBeginning()
    txt = ['<MRML  version="Slicer4.4.0" userTags="">',...
     '<Crosshair id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="JumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"></Crosshair>\r\n',...
     '<Selection id="vtkMRMLSelectionNodeSingleton" name="Selection" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" frequencyUnitNodeRef="vtkMRMLUnitNodeApplicationFrequency" intensityUnitNodeRef="vtkMRMLUnitNodeApplicationIntensity" lengthUnitNodeRef="vtkMRMLUnitNodeApplicationLength" timeUnitNodeRef="vtkMRMLUnitNodeApplicationTime" velocityUnitNodeRef="vtkMRMLUnitNodeApplicationVelocity" references="unit/frequency:vtkMRMLUnitNodeApplicationFrequency;unit/intensity:vtkMRMLUnitNodeApplicationIntensity;unit/length:vtkMRMLUnitNodeApplicationLength;unit/time:vtkMRMLUnitNodeApplicationTime;unit/velocity:vtkMRMLUnitNodeApplicationVelocity;" activeVolumeID="vtkMRMLScalarVolumeNode2" secondaryVolumeID="NULL" activeLabelVolumeID="NULL" activeFiducialListID="NULL" activePlaceNodeID="NULL" activePlaceNodeClassName="NULL" activeROIListID="NULL" activeCameraID="NULL" activeTableID="NULL" activeViewID="NULL" activeLayoutID="NULL" ></Selection>\r\n',...
     '<Interaction id="vtkMRMLInteractionNodeSingleton" name="Interaction" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" currentInteractionMode="ViewTransform" placeModePersistence="false" lastInteractionMode="ViewTransform" ></Interaction>\r\n',...
     '<View id="vtkMRMLViewNode1" name="View1" hideFromEditors="false" selectable="true" selected="false" singletonTag="1" attributes="MappedInLayout:1" layoutLabel="1" layoutName="1" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="200" letterSize="0.05" boxVisible="true" fiducialsVisible="true" fiducialLabelsVisible="true" axisLabelsVisible="true" axisLabelsCameraDependent="true" animationMode="Off" viewAxisMode="LookFrom" spinDegrees="2" spinMs="5" spinDirection="YawLeft" rotateDegrees="5" rockLength="200" rockCount="0" stereoType="NoStereo" renderMode="Perspective" useDepthPeeling="0" ></View>\r\n',...
     '<Slice id="vtkMRMLSliceNodeRed" name="Red" hideFromEditors="false" selectable="true" selected="false" singletonTag="Red" attributes="MappedInLayout:1" layoutLabel="R" layoutName="Red" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="318.002 259.289 0.5" dimensions="1002 817 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="318.002 259.289 0.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 1 0 -1.77324e-006 0 0 1 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.952941 0.290196 0.2" orientation="Axial" orientationReference="Axial" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n',...
     '<Slice id="vtkMRMLSliceNodeYellow" name="Yellow" hideFromEditors="false" selectable="true" selected="false" singletonTag="Yellow" attributes="MappedInLayout:1" layoutLabel="Y" layoutName="Yellow" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="0 0 1 0.218502 -1 0 0 -1.77324e-006 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.929412 0.835294 0.298039" orientation="Sagittal" orientationReference="Sagittal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n',...
     '<Slice id="vtkMRMLSliceNodeGreen" name="Green" hideFromEditors="false" selectable="true" selected="false" singletonTag="Green" attributes="MappedInLayout:1" layoutLabel="G" layoutName="Green" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 0 1 0.218498 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.431373 0.690196 0.294118" orientation="Coronal" orientationReference="Coronal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1" ></Slice>\r\n',...
     '<Layout id="vtkMRMLLayoutNodevtkMRMLLayoutNode" name="Layout" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLLayoutNode" currentViewArrangement="3" guiPanelVisibility="1" bottomPanelVisibility ="1" guiPanelLR="0" collapseSliceControllers="0"',...
     ' numberOfCompareViewRows="1" numberOfCompareViewColumns="1" numberOfLightboxRows="6" numberOfLightboxColumns="6" mainPanelSize="400" secondaryPanelSize="400" ></Layout>\r\n',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeRed" name="SliceComposite" hideFromEditors="true" selectable="true" selected="false" singletonTag="Red" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Red" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeYellow" name="SliceComposite_1" hideFromEditors="true" selectable="true" selected="false" singletonTag="Yellow" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Yellow" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n',...
     '<SliceComposite id="vtkMRMLSliceCompositeNodeGreen" name="SliceComposite_2" hideFromEditors="true" selectable="true" selected="false" singletonTag="Green" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Green" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1" ></SliceComposite>\r\n',...
     '<Camera id="vtkMRMLCameraNode1" name="Default Scene Camera" hideFromEditors="false" selectable="true" selected="false" userTags="" position="-174.789 463.886 65.26" focalPoint="0 0 0" viewUp="0.0301267 -0.128106 0.991303" parallelProjection="false" parallelScale="1" activetag="vtkMRMLViewNode1" appliedTransform="1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1" ></Camera>\r\n',...
     '<SubjectHierarchy id="vtkMRMLSubjectHierarchyNode1" name="SubjectHierarchy" hideFromEditors="false" selectable="true" selected="false" attributes="SubjectHierarchyVersion:2" >\r\n',...
       '<SubjectHierarchyItem id="3" name="Scene" parent="0" type="" expanded="true" attributes="Level^Scene|">\r\n',...
         '<SubjectHierarchyItem id="9" dataNode="vtkMRMLScalarVolumeNode1" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>\r\n',...
         '<SubjectHierarchyItem id="12" dataNode="vtkMRMLScalarVolumeNode2" parent="3" type="Volumes" expanded="true"></SubjectHierarchyItem>\r\n',...
         '<SubjectHierarchyItem id="13" dataNode="vtkMRMLSceneViewNode1" parent="3" type="SceneViews" expanded="true"></SubjectHierarchyItem>\r\n',...
         '<SubjectHierarchyItem id="15" expanded="true" type="Markups" parent="3" dataNode="vtkMRMLMarkupsFiducialNode1"></SubjectHierarchyItem></SubjectHierarchyItem></SubjectHierarchy>\r\n',...
     '<ClipModels id="vtkMRMLClipModelsNodevtkMRMLClipModelsNode" name="ClipModels" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLClipModelsNode" clipType="0" redSliceClipState="0" yellowSliceClipState="0" greenSliceClipState="0" ></ClipModels>\r\n',...
     '<ScriptedModule id="vtkMRMLScriptedModuleNodeDataProbe" name="ScriptedModule" hideFromEditors="true" selectable="true" selected="false" singletonTag="DataProbe" ModuleName ="DataProbe" ></ScriptedModule>\r\n'...
     ];
end

function txt = GetFiducialEnding(fiducialFileName)
    [~,name_no_ext_no_path,~]=fileparts(fiducialFileName);
    txt = ['<Crosshair id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="JumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"  ></Crosshair>',...
    '<Selection id="vtkMRMLSelectionNodeSingleton" name="Selection" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" frequencyUnitNodeRef="vtkMRMLUnitNodeApplicationFrequency" intensityUnitNodeRef="vtkMRMLUnitNodeApplicationIntensity" lengthUnitNodeRef="vtkMRMLUnitNodeApplicationLength" timeUnitNodeRef="vtkMRMLUnitNodeApplicationTime" velocityUnitNodeRef="vtkMRMLUnitNodeApplicationVelocity" references="unit/frequency:vtkMRMLUnitNodeApplicationFrequency;unit/intensity:vtkMRMLUnitNodeApplicationIntensity;unit/length:vtkMRMLUnitNodeApplicationLength;unit/time:vtkMRMLUnitNodeApplicationTime;unit/velocity:vtkMRMLUnitNodeApplicationVelocity;" activeVolumeID="vtkMRMLScalarVolumeNode2" secondaryVolumeID="vtkMRMLScalarVolumeNode2" activeLabelVolumeID="NULL" activeFiducialListID="NULL" activePlaceNodeID="NULL" activePlaceNodeClassName="NULL" activeROIListID="NULL" activeCameraID="NULL" activeTableID="NULL" activeViewID="NULL" activeLayoutID="NULL"  ></Selection>',...
    '<Interaction id="vtkMRMLInteractionNodeSingleton" name="Interaction" hideFromEditors="true" selectable="true" selected="false" singletonTag="Singleton" currentInteractionMode="ViewTransform" placeModePersistence="false" lastInteractionMode="ViewTransform"  ></Interaction>',...
    '<View id="vtkMRMLViewNode1" name="View1" hideFromEditors="false" selectable="true" selected="false" singletonTag="1" attributes="MappedInLayout:1" layoutLabel="1" layoutName="1" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="200" letterSize="0.05" boxVisible="true" fiducialsVisible="true" fiducialLabelsVisible="true" axisLabelsVisible="true" axisLabelsCameraDependent="true" animationMode="Off" viewAxisMode="LookFrom" spinDegrees="2" spinMs="5" spinDirection="YawLeft" rotateDegrees="5" rockLength="200" rockCount="0" stereoType="NoStereo" renderMode="Perspective" useDepthPeeling="0"  ></View>',...
    '<Slice id="vtkMRMLSliceNodeRed" name="Red" hideFromEditors="false" selectable="true" selected="false" singletonTag="Red" attributes="MappedInLayout:1" layoutLabel="R" layoutName="Red" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="318.002 259.289 0.5" dimensions="1002 817 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="318.002 259.289 0.5" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 1 0 -1.77324e-006 0 0 1 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.952941 0.290196 0.2" orientation="Axial" orientationReference="Axial" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Slice id="vtkMRMLSliceNodeYellow" name="Yellow" hideFromEditors="false" selectable="true" selected="false" singletonTag="Yellow" attributes="MappedInLayout:1" layoutLabel="Y" layoutName="Yellow" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="0 0 1 0.218502 -1 0 0 -1.77324e-006 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.929412 0.835294 0.298039" orientation="Sagittal" orientationReference="Sagittal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Slice id="vtkMRMLSliceNodeGreen" name="Green" hideFromEditors="false" selectable="true" selected="false" singletonTag="Green" attributes="MappedInLayout:1" layoutLabel="G" layoutName="Green" active="false" visibility="true" backgroundColor="0 0 0" backgroundColor2="0 0 0" orientationMarkerType="none" orientationMarkerSize="medium" rulerType="none" AxisLabels="L;R;P;A;I;S" fieldOfView="264.021 215.011 0.437" dimensions="1002 816 1" xyzOrigin="0 0 0" sliceResolutionMode="1" uvwExtents="264.021 215.011 0.437" uvwDimensions="256 256 1" uvwOrigin="0 0 0" activeSlice="0" layoutGridRows="1" layoutGridColumns="1" sliceToRAS="-1 0 0 1.77324e-006 0 0 1 0.218498 0 1 0 0 0 0 0 1" orientationMatrixAxial="-1 0 0 0 1 0 0 0 1" orientationMatrixSagittal="0 0 1 -1 0 0 0 1 0" orientationMatrixCoronal="-1 0 0 0 0 1 0 1 0" layoutColor="0.431373 0.690196 0.294118" orientation="Coronal" orientationReference="Coronal" jumpMode="1" sliceVisibility="true" widgetVisibility="false" useLabelOutline="false" sliceSpacingMode="0" prescribedSliceSpacing="1 1 1"  ></Slice>',...
    '<Layout id="vtkMRMLLayoutNodevtkMRMLLayoutNode" name="Layout" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLLayoutNode" currentViewArrangement="3" guiPanelVisibility="1" bottomPanelVisibility ="1" guiPanelLR="0" collapseSliceControllers="0"',...
    ' numberOfCompareViewRows="1" numberOfCompareViewColumns="1" numberOfLightboxRows="6" numberOfLightboxColumns="6" mainPanelSize="400" secondaryPanelSize="400"  ></Layout>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeRed" name="SliceComposite" hideFromEditors="true" selectable="true" selected="false" singletonTag="Red" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Red" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeYellow" name="SliceComposite_1" hideFromEditors="true" selectable="true" selected="false" singletonTag="Yellow" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Yellow" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<SliceComposite id="vtkMRMLSliceCompositeNodeGreen" name="SliceComposite_2" hideFromEditors="true" selectable="true" selected="false" singletonTag="Green" backgroundVolumeID="vtkMRMLScalarVolumeNode2" foregroundVolumeID="vtkMRMLScalarVolumeNode1" labelVolumeID="" compositing="0" foregroundOpacity="0.5" labelOpacity="1" linkedControl="1" fiducialVisibility="1" fiducialLabelVisibility="1" sliceIntersectionVisibility="0" layoutName="Green" annotationSpace="IJKAndRAS" annotationMode="All" doPropagateVolumeSelection="1"  ></SliceComposite>',...
    '<Camera id="vtkMRMLCameraNode1" name="Default Scene Camera" hideFromEditors="false" selectable="true" selected="false" userTags="" position="-174.789 463.886 65.26" focalPoint="0 0 0" viewUp="0.0301267 -0.128106 0.991303" parallelProjection="false" parallelScale="1" activetag="vtkMRMLViewNode1" appliedTransform="1 0 0 0  0 1 0 0  0 0 1 0  0 0 0 1"  ></Camera>',...
    '<ClipModels id="vtkMRMLClipModelsNodevtkMRMLClipModelsNode" name="ClipModels" hideFromEditors="true" selectable="true" selected="false" singletonTag="vtkMRMLClipModelsNode" clipType="0" redSliceClipState="0" yellowSliceClipState="0" greenSliceClipState="0"  ></ClipModels>',...
    '<ScriptedModule id="vtkMRMLScriptedModuleNodeDataProbe" name="ScriptedModule" hideFromEditors="true" selectable="true" selected="false" singletonTag="DataProbe" ModuleName ="DataProbe"  ></ScriptedModule>\r\n',...
    '<SceneViewStorage id="vtkMRMLSceneViewStorageNode1" name="SceneViewStorage" hideFromEditors="true" selectable="true" selected="false" fileName="MasterSceneView.png" useCompression="1" defaultWriteFileExtension="png" readState="4" writeState="4" ></SceneViewStorage>\r\n',...
    '<MarkupsFiducialStorage selected="false" selectable="true" hideFromEditors="true" name="MarkupsFiducialStorage" id="vtkMRMLMarkupsFiducialStorageNode1" writeState="0" readState="0" defaultWriteFileExtension="fcsv" useCompression="1" fileName="',fiducialFileName,'" coordinateSystem="RAS"></MarkupsFiducialStorage>\r\n',...
    '<MarkupsFiducial userTags="" selected="false" selectable="true" hideFromEditors="false" name="' name_no_ext_no_path '" id="vtkMRMLMarkupsFiducialNode1" references="display:vtkMRMLMarkupsFiducialDisplayNode1;storage:vtkMRMLMarkupsFiducialStorageNode1;" storageNodeRef="vtkMRMLMarkupsFiducialStorageNode1" displayNodeRef="vtkMRMLMarkupsFiducialDisplayNode1" markupLabelFormat="%%N-%%d" locked="0"></MarkupsFiducial>\r\n',...
    '<MarkupsFiducialDisplay   id="vtkMRMLMarkupsFiducialDisplayNode1" name="MarkupsFiducialDisplay" hideFromEditors="true" selectable="true" selected="false" color="0.4 1 0" edgeColor="0 0 0" selectedColor="1 0.500008 0.500008" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" visibility2D="true" visibility3D="true" edgeVisibility="false" clipping="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="false" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" activeAttributeLocation="point" viewNodeRef="" folderDisplayOverrideAllowed="true" propertiesLabelVisibility="false" pointLabelsVisibility="true" textScale="3" glyphScale="1" glyphSize="5" useGlyphScale="true" glyphType="Sphere3D" snapMode="toVisibleSurface" sliceProjection="false" sliceProjectionUseFiducialColor="true" sliceProjectionOutlinedBehindSlicePlane="false" sliceProjectionColor="1 1 1" sliceProjectionOpacity="0.6" curveLineSizeMode="UseLineThickness" lineThickness="0.2" lineDiameter="1" lineColorFadingStart="1" lineColorFadingEnd="10" lineColorFadingSaturation="1" lineColorFadingHueOffset="0" handlesInteractive="false" translationHandleVisibility="true" rotationHandleVisibility="true" scaleHandleVisibility="true" fillVisibility="true" outlineVisibility="true" fillOpacity="0.5" outlineOpacity="1" occludedVisibility="false" occludedOpacity="0.3" textProperty="font-family:Arial;font-size:5px;font-style:normal;font-weight:normal;color:rgba(255,255,255,1);background-color:rgba(0,0,0,0);border-width:1px;border-color:rgba(255,255,255,0.0);text-shadow:1px 1px 0px rgba(0,0,0,0.0);" activeColor="0.4 1 0" ></MarkupsFiducialDisplay>\r\n'...
    '</MRML>'
    ];
end
