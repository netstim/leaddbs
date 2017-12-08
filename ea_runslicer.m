%% Function to launch slicer and load *.nii files
%  Last Revision: 7/12/2017 
%  Thushara Perera (c) 2017 Bionics Institute
%  Input:
%   - lead dbs options struct
%   - task integer ID to tell function which files to open 
%     (see inline comments for more detail)
%  Output:
%   - relevant slicer scene (mrml) file will be saved in patient folder
%   - temporary python script file will be written to lead_dbs root
%
%  Things to note: user must set the path of the slicer executable. It would
%  be nice if there was an easy method to determine the path automatically.
%  
%% 
function ea_runslicer(options, task)
    options.prefs = ea_prefs('');
    slicer_path = options.prefs.slicer.dir;
    [pth,fn,ext]=fileparts(slicer_path);
    if strcmp(ext,'.app') && ismac
        slicer_path=fullfile(slicer_path,'Contents','MacOS','Slicer');
    end
    lead_path = options.earoot;
    slicer_mrml = 'none';
    
    if (~exist(slicer_path, 'file'))
        warning('Path to Slicer executable not found');
        return;
    end
    
    if (isempty(options.uipatdirs))
        warning('No patient selected');
        return;
    else
        patient_path = options.uipatdirs{1};
    end
    
    switch task
        case 1 % show unnormalised
            nfiles = 0;
            slicer_mrml = 'Slicer_unnormalised.mrml';
            allfiles = {options.prefs.tranii_unnormalized
                options.prefs.sagnii_unnormalized
                options.prefs.cornii_unnormalized
                options.prefs.rawctnii_unnormalized
                options.prefs.prenii_unnormalized
                options.prefs.prenii_unnormalized_t1
                options.prefs.prenii_unnormalized_pd};
            
            for i=1:length(allfiles)
                if (exist([patient_path, filesep, allfiles{i}], 'file') == 2)
                    nfiles = nfiles + 1;
                    filenames{nfiles} = allfiles{i};
                    filepaths{nfiles} = [patient_path, filesep, allfiles{i}];
                end
            end
            
        case 2 % show data after co-registration
            nfiles = 0;
            slicer_mrml = 'Slicer_coregistered.mrml';
            allfiles = {options.prefs.tranii_unnormalized
                options.prefs.sagnii_unnormalized
                options.prefs.cornii_unnormalized
                options.prefs.ctnii_coregistered % check others
                options.prefs.prenii_unnormalized
                options.prefs.prenii_unnormalized_t1
                options.prefs.prenii_unnormalized_pd
                };
            
            for i=1:length(allfiles)
                if (exist([patient_path, filesep, allfiles{i}], 'file') == 2)
                    nfiles = nfiles + 1;
                    filenames{nfiles} = allfiles{i};
                    filepaths{nfiles} = [patient_path, filesep, allfiles{i}];
                end
            end
            
        case 3 % show data after normalisation
            nfiles = 0;
            slicer_mrml = 'Slicer_normalised.mrml';
            allfiles = {
                't2.nii'
                'glanat_t2.nii'   % could not find pref for this in options.prefs.gctnii
                options.prefs.gprenii
                options.prefs.gctnii
                options.prefs.gtranii
                options.prefs.gcornii
                options.prefs.gsagnii
                };
            
            if (exist(allfiles{1}, 'file') == 2) %special case for template
                nfiles = nfiles + 1;
                filenames{nfiles} = allfiles{1};
                 % is there a better way to get the MNI template?
                filepaths{nfiles} = [lead_path, 'templates', filesep, 'space', filesep, 'MNI_ICBM_2009b_NLIN_ASYM', filesep, allfiles{1}];
            end
            
            for i=2:length(allfiles)
                if (exist([patient_path, filesep, allfiles{i}], 'file') == 2)
                    nfiles = nfiles + 1;
                filenames{nfiles} = allfiles{i};
                    filepaths{nfiles} = [patient_path, filesep, allfiles{i}];
                end
            end
            
        otherwise
            warning('Task ID not recognised');
            return;
    end

    script_path = [lead_path, 'slicer.py'];
    scene_path = [patient_path, filesep, slicer_mrml];
    if (nfiles < 1)
        warning('Need at least one volume image to load into slicer');
        return;
    end
    
    fid = fopen(scene_path, 'w');
    fprintf(fid, [GetBeginning(), '\r\n\r\n']);
    for i=1:nfiles
        fprintf(fid, [GetFileXML(i, filepaths{i}, filenames{i}), '\r\n\r\n']);
    end
    fprintf(fid, GetEnding());
    fclose(fid);
    
    fid = fopen(script_path,'w'); % write temp python script to load volumes
    fprintf(fid, ['slicer.util.loadScene("', strrep(scene_path, '\', '\\'), '")\r\n']); 
    fclose(fid);
    disp('Loading up 3D Slicer...');
    system(['"', slicer_path, '" --no-splash --python-script "', script_path, '" &']); 
    % the trailing '&' returns control back to matlab without waiting for slicer to close
end

function txt = GetFileXML(index, filepath, name)

    path = strrep(filepath, '\', '\\'); %avoid escape char errors

    vdisplay = ['<VolumeDisplay id="vtkMRMLScalarVolumeDisplayNode', num2str(index), '" ',...
        'name="VolumeDisplay" hideFromEditors="true" selectable="true" selected="false" color="0.5 0.5 0.5" edgeColor="0 0 0" selectedColor="1 0 0" selectedAmbient="0.4" ambient="0" diffuse="1" selectedSpecular="0.5" specular="0" power="1" opacity="1" sliceIntersectionOpacity="1" pointSize="1" lineWidth="1" representation="2" lighting="true" interpolation="1" shading="true" visibility="true" edgeVisibility="false" clipping="false" sliceIntersectionVisibility="false" sliceIntersectionThickness="1" frontfaceCulling="false" backfaceCulling="true" scalarVisibility="false" vectorVisibility="false" tensorVisibility="false" interpolateTexture="false" scalarRangeFlag="UseData" scalarRange="0 100" colorNodeID="vtkMRMLColorTableNodeGrey"  window="100" level="50" upperThreshold="32767" lowerThreshold="-32768" interpolate="1" autoWindowLevel="0" applyThreshold="0" autoThreshold="0" ></VolumeDisplay>'];


    vstorage = ['<VolumeArchetypeStorage ',...
        'id="vtkMRMLVolumeArchetypeStorageNode', num2str(index),'" ',...
        'name="VolumeArchetypeStorage_', num2str(index),'" ',...
        'hideFromEditors="true" selectable="true" selected="false" ',...
        'fileName="',path,'" ',...
        'useCompression="1" defaultWriteFileExtension="nrrd" readState="0" writeState="0" centerImage="0" UseOrientationFromFile="1" ></VolumeArchetypeStorage>'];


    volume = ['<Volume id="vtkMRMLScalarVolumeNode', num2str(index),'" ',...
        'name="',name(1:end-4),'" ',... % remove extension .nii
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
    txt = ['<SceneView id="vtkMRMLSceneViewNode1" name="Master Scene View" hideFromEditors="false" selectable="true" selected="false" storageNodeRef="vtkMRMLSceneViewStorageNode1" references="storage:vtkMRMLSceneViewStorageNode1;" userTags="" screenshotType="4" sceneViewDescription="Scene at MRML file save point" >  <Crosshair',...
    ' id="vtkMRMLCrosshairNodedefault" name="Crosshair" hideFromEditors="true" selectable="true" selected="false" singletonTag="default" crosshairMode="NoCrosshair" navigation="false" crosshairBehavior="JumpSlice" crosshairThickness="Fine" crosshairRAS="0 0 0"  ></Crosshair>',...
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
    '</SceneView>',...
    '<SceneViewStorage id="vtkMRMLSceneViewStorageNode1" name="SceneViewStorage" hideFromEditors="true" selectable="true" selected="false" fileName="Master Scene View.png" useCompression="1" defaultWriteFileExtension="png" readState="0" writeState="4" ></SceneViewStorage>',...
    '</MRML>'
    ];
end