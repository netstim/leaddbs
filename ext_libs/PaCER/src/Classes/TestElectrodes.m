%% TestElectrodes - Class representing a Ben's Gun configuration of Test Electrodes (up to rotation!)
% note that the precise rotation depends on the kinematics of the
% sterotactic frame/the micro drive device
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
classdef TestElectrodes < Trajectory & plotable3D & plotableImage & configurable & id
    properties (SetAccess = protected, GetAccess=public, SetObservable = true) 
       originalEntryPoint = NaN(3,1);
       originalTargetPoint = NaN(3,1);
       marker = 'x';       %default circle
       color = [0.2 1 0.4];    %default 
       diameter = 4;       % groesse
       lineStyle = '--';      
       show2DProjection = 1;  % variable {0,1} die angibt, ob 2D plots angezeigt werden oder nicht
       %INHERITED
       %plotHandleHashMap (from plotableImage)
       trajectoryChangedListener = event.listener.empty();
    end
    
    properties(Access=public, SetObservable, AbortSet)
        obliqueViewRadius = 80; 
        currentDepthValue = 0; % 0= we ar at target; -x means we are x mm above the target; +x means the depth is x + target in mm
    end
    
    properties (SetAccess = public, GetAccess=public, SetObservable = true) 
        intersections = zeros(2,5) ; %TODO: [tg1, tg2] for all electrodes %center, medial, lateral,  posterior, anterior. Set  NaN for unused Electrodes. (Replace by Dep. Prop. that filtes NaN)
        ELECTRODE_DISTANCE = 2.0; %[mm], distance from posterior, anterior, lateral, medial to center electrode
        ELECTRODE_RECORDED_INTERSECTION_COLOR = 'r';
        ELECTRODE_INTERSECTION_COLOR = 'g';
        ELECTRODE_COLOR = [0 1 1]; %lightblue
        ELECTRODE_PLOT_EXTENSION_MM = 15; %[mm], define how long the plottet electrode lines should be extended beyond the target/entry points;
        DISPLAY_TARGET_POINTS = false;
    end
    
    properties (SetAccess = protected, GetAccess = public, ... %dependent props.
                Dependent = true) 
       electrodesEntryPoints = {}; %center, right, left,  posterior, anterior.
       electrodesTargetPoints = {}; %center, right, left,  posterior, anterior.
    end
    properties (Access = private)
       electrodesOffsets;
       electrodesOffsetsDirection;
    end
    properties(SetAccess = private, GetAccess = public, Dependent = true)
        % to control ObliqueView
        trajectoryLength = 0;
    end
    events
       SliderValueChanged; %OBSOLET
    end
    methods
        function this = TestElectrodes(varargin)
            if(nargin == 2)
                entryPoint = varargin{1};
                targetPoint = varargin{2};
                
                if(sum(size(entryPoint) ~= [3 1]) ||sum(size(targetPoint) ~= [3 1]))
                    warning('TestElectrodes: entryPoint and TargetPoint MUST be 3x1 vectors!'); %#ok<WNTAG>
                    if(size(entryPoint) ~= [3 1]) %#ok<BDSCA>
                        entryPoint = entryPoint';
                    end
                    if(size(targetPoint) ~= [3 1]) %#ok<BDSCA>
                        targetPoint = targetPoint';
                    end
                end
                this.entryPoint = entryPoint;
                this.targetPoint = targetPoint;
                this.originalEntryPoint = this.entryPoint;
                this.originalTargetPoint = this.targetPoint;
                this.currentDepthValue = -this.trajectoryLength;
            end
            addlistener(this, 'ELECTRODE_COLOR', 'PostSet', @this.updatePlot3D);
            this.trajectoryChangedListener = this.addlistener('trajectoryChanged', @this.onUpdatePlots);
        end
        
        function electrodesEntryPoints = get.electrodesEntryPoints(this)
            centralpoint = repmat(this.entryPoint, 1, 5); % x y z for 5 electrodes
            offsets = this.getElectrodesOffsets();
            electrodesEntryPoints = centralpoint + offsets;
        end
 
        function electrodesTargetPoints = get.electrodesTargetPoints(this)
            centralpoint = repmat(this.targetPoint, 1, 5); % x y z for 5 electrodes
            offsets = this.getElectrodesOffsets();
            electrodesTargetPoints = centralpoint + offsets;
                                 %   electrodesTargetPoints = mat2cell(centralpoint + offsets, ...
                                %     3,ones(length(centralpoint),1));
        end
        
        function trajectoryLength = get.trajectoryLength(this)
            trajectoryLength = this.entryPoint-this.targetPoint;
            trajectoryLength = ceil(sqrt(dot(trajectoryLength, trajectoryLength)));
        end
        
        function obliqueViewRadius = get.obliqueViewRadius(this)
            obliqueViewRadius = this.obliqueViewRadius;
        end
        
        function set.obliqueViewRadius(this, value)
            if(value < 0)
                disp(strcat('Testelectrodes:set.obliqueViewRadius:Radius have to be greater than 0'));
            else
                this.obliqueViewRadius = value;
            end;
        end        
        
        function pointCloud = getIntersectionPointCloud(this)
            % 3d-point = Target + rg * tg1          
            pointCloudTg1 = this.electrodesTargetPoints + repmat(this.direction,1,5) .* repmat(this.intersections(1,:),3,1);
            
            % 3d-point = Target + rg * tg2
            pointCloudTg2 = this.electrodesTargetPoints + repmat(this.direction,1,5) .* repmat(this.intersections(2,:),3,1);
            
            pointCloud = [pointCloudTg1 pointCloudTg2]; % TODO Remove Electrodes with no Intersection 
            
        end
        
        function disableTrajectoryChangedEvent(this)
            this.trajectoryChangedListener.delete();
        end
        
        function enableTrajectoryChangedEvent(this)
            this.trajectoryChangedListener = this.addlistener('trajectoryChanged', @this.onUpdatePlots);
        end
    end
    
    methods (Access = private)
        function [offsets, B] = getElectrodesOffsets(this)
            if(~isempty(this.electrodesOffsets) && ~isnan(this.electrodesOffsets(1,1)) &&isequal(this.direction, this.electrodesOffsetsDirection))
                offsets = this.electrodesOffsets;
            else
                config = this.ELECTRODE_DISTANCE * ...
                    [0 0 0; % 5-dice configuration, FIXME: should be ortogonal 2mm to center electrode
                    0 1 0; % note that config is transposed!
                    0 -1 0;
                    0 0 -1;
                    0 0 1];
                
                e1 = [1; 0; 0]; % <-> l-r TODO: replace this by direction vecotrs from frame
                e2 = [0; 1; 0]; % <-> a-p
                
                b1 = this.direction; %MUST be normalized!
                if(all(eq(abs(round(b1)), [1;0;0])))
                    b2  = [0; 1; 0];
                    b3 = [0; 0; 1];
                else
                    b2 = e1 - dot(b1,e1) * b1;
                    b2 = b2 / norm(b2);
                    
                    b3 = e2 - dot(b1, e2) * b1 - dot(b2,e2) * b2;
                    b3 = b3 / norm(b3);
                end
                
                
                
                B = [b1,b2,b3]; %new BASE (b2,b3 forming a plane ortogonal to our electrode)
                offsets = B * config';
                this.electrodesOffsets = offsets;
                this.electrodesOffsetsDirection = this.direction;
            end
        end
        
        function [base,directionVec] = getElectrodePoints(this)
            base = this.electrodesTargetPoints;
            directionVec = (this.electrodesEntryPoints-base);
        end%fcn-getElectrodePoints
        
    end
    
    methods % graphical methods
        function graphicsHandle = initPlot3D(this, parentAxes)
            if(this.plotHandleHashMap3D.isKey(parentAxes))
                warning('TestElectrodes:initPlot3DcalledMoreThanOnce', 'initPlot3D was called more than once for the same axes and the same object!');
                graphicsHandle = this.plotHandleHashMap3D(parentAxes);
                return;
            end
            % set(0,'DefaultLineLineSmoothing','on'); %enable line smoothing for all objects that support it
            
            % create a group object and group all plots to this "parent" handle
            graphicsHandle = hggroup('Parent', parentAxes);
            
            d = this.electrodesTargetPoints;
            e = this.electrodesEntryPoints;
            tg = this.getIntersectionPointCloud;
            
            set(parentAxes, 'NextPlot', 'add');
            if(this.DISPLAY_TARGET_POINTS)
                graphicsHandles.targetPoints = plot3(d(1,:),d(2,:),d(3,:), 'b*', 'Parent', graphicsHandle,'Clipping','on');
                set(graphicsHandles.targetPoints, 'Parent', graphicsHandle);
            end
            graphicsHandles.entryPoints = plot3(e(1,:),e(2,:),e(3,:), 'g*', 'Parent', graphicsHandle,'Clipping','on');
            set(graphicsHandles.entryPoints, 'Parent', graphicsHandle);
            graphicsHandles.intersectionPoints = plot3(tg(1,:),tg(2,:),tg(3,:), 'r*', 'Parent', graphicsHandle,'Clipping','on');
            set(graphicsHandles.intersectionPoints, 'Parent', graphicsHandle);
            
            % plot trajectories XX mm beyoned the target point and before the entry point
            dPlusXX = d + repmat(this.direction,1,5) * this.ELECTRODE_PLOT_EXTENSION_MM;
            eMinusXX = e - repmat(this.direction,1,5) * this.ELECTRODE_PLOT_EXTENSION_MM;
            
            graphicsHandles.electrodes = ...
                line([eMinusXX(1,:); dPlusXX(1,:)], ...
                [eMinusXX(2,:); dPlusXX(2,:)], ...
                [eMinusXX(3,:); dPlusXX(3,:)], ...
                'LineWidth', 1.0,...
                'Color', this.ELECTRODE_COLOR,...
                'Parent', graphicsHandle);
            set(graphicsHandles.electrodes , 'Parent', graphicsHandle);
            
            % mark lateral (left)
            graphicsHandles.electrodes = ...
                line([eMinusXX(1,3); dPlusXX(1,3)], ...
                [eMinusXX(2,3); dPlusXX(2,3)], ...
                [eMinusXX(3,3); dPlusXX(3,3)], ...
                'LineWidth', 1.0,...
                'Color', this.ELECTRODE_COLOR,...
                'Parent', graphicsHandle);
            set(graphicsHandles.electrodes , 'Parent', graphicsHandle);
            
            graphicsHandles.recordedIntersections = ...
                line([tg(1,1:5); tg(1,6:10)], ...
                [tg(2,1:5); tg(2,6:10)], ...
                [tg(3,1:5); tg(3,6:10)], ...
                'LineWidth', 3.0,...
                'Color', this.ELECTRODE_RECORDED_INTERSECTION_COLOR ,...
                'Parent', graphicsHandle);
            set(graphicsHandles.recordedIntersections , 'Parent', graphicsHandle);
            
            this.plotHandleHashMap3D(double(parentAxes)) = double(graphicsHandle);
            
        end%fcn-initPlot3D
        
%TODO: dragging points implementieren        
        function graphicsHandle = initPlotAxial(this, ax, CurrentPositionObject)
            % create a group object and group all plots to this "parent" handle
            graphicsHandle = hggroup('Parent', ax);
%             keyVal = ['axial',this.String,ax];
            keyVal = ax;
%             [base, direcVec] = this.getElectrodePoints;
%             curPosInMM = repmat(CurrentPositionObject.currentImagePosInMmA,1,5);
%             factor = repmat(((curPosInMM - base(3,:))./direcVec(3,:)),3,1);
            curPosInMM = CurrentPositionObject.currentImagePosInMmA;
            base = this.targetPoint;
            direcVec = this.entryPoint - base;
            factor = ((curPosInMM - base(3))./direcVec(3));
            projectionPoint = base + factor.*direcVec;
            %if abfrage wenn projection < targetPoint, dann p = t oder
            %keine projektion
            pointStart = [this.entryPoint(1) this.targetPoint(1)];
            pointEnd = [this.entryPoint(2) this.targetPoint(2)];

            projectionPointStart = [this.entryPoint(1) projectionPoint(1)];
            projectionPointEnd = [this.entryPoint(2) projectionPoint(2)];
            
            posLineColor = [0 1 1];
            
            %check wheather we are within entry-target or >taget
            if(projectionPoint(2) < min(this.targetPoint(2),this.entryPoint(2)))%note it starts from left-bottom-corner
                projectionPointStart = [this.targetPoint(1) projectionPoint(1)];
                projectionPointEnd = [this.targetPoint(2) projectionPoint(2)];
                posLineColor = [1 0 0];
            elseif(projectionPoint(2) > max(this.entryPoint(2),this.targetPoint(2)))
                projectionPointStart = [this.entryPoint(1) projectionPoint(1)];
                projectionPointEnd = [this.entryPoint(2) projectionPoint(2)];
                posLineColor = [1 0 0];
            end%if
            %start to plot what needs to be plotted
            handle.CurrentPositionObject = CurrentPositionObject;
%             this.plotHandleHashMap(1.0) = CurrentPositionObject;
            this.onOrthogonalProjectionPlot(keyVal, graphicsHandle, pointStart, pointEnd, projectionPointStart, projectionPointEnd, posLineColor, handle);
        
        end%fcn-initPlotAxial
        
        function graphicsHandle = initPlotCoronal(this, ax, CurrentPositionObject)
            % create a group object and group all plots to this "parent" handle
            graphicsHandle = hggroup('Parent', ax);
%             keyVal = ['coronal',this.String,ax];
            keyVal = ax;
%             [base, direcVec] = this.getElectrodePoints;
%             curPosInMM = repmat(CurrentPositionObject.currentImagePosInMmC,1,5);
%             factor = repmat(((curPosInMM - base(2,:))./direcVec(2,:)),3,1);
%             projectionPoint = base + factor.*direcVec;
            curPosInMM = CurrentPositionObject.currentImagePosInMmC;
            base = this.targetPoint;
            direcVec = this.entryPoint - base;
            factor = ((curPosInMM - base(2))./direcVec(2));
            projectionPoint = base + factor.*direcVec;
            
            pointStart = [this.entryPoint(1) this.targetPoint(1)];
            pointEnd = [this.entryPoint(3) this.targetPoint(3)];
            
            projectionPointStart = [this.entryPoint(1) projectionPoint(1)];
            projectionPointEnd = [this.entryPoint(3) projectionPoint(3)];
            
            posLineColor = [0 1 1];
            
            %check wheather we are within entry-target or >taget
            if(projectionPoint(3) < min(this.targetPoint(3),this.entryPoint(3)))%note it starts from left-bottom-corner
                projectionPointStart = [this.targetPoint(1) projectionPoint(1)];
                projectionPointEnd = [this.targetPoint(3) projectionPoint(3)];
                posLineColor = [1 0 0];
            elseif(projectionPoint(3) > max(this.entryPoint(3),this.targetPoint(3)))
                projectionPointStart = [this.entryPoint(1) projectionPoint(1)];
                projectionPointEnd = [this.entryPoint(3) projectionPoint(3)];
                posLineColor = [1 0 0];
            end%if
            %start to plot what needs to be plotted
            handle.CurrentPositionObject = CurrentPositionObject;
%             this.plotHandleHashMap(1.0) = CurrentPositionObject;
            this.onOrthogonalProjectionPlot(keyVal, graphicsHandle, pointStart, pointEnd, projectionPointStart, projectionPointEnd, posLineColor, handle);
        
        end%fcn-initPlotCoronal
        
        function graphicsHandle = initPlotSagital(this,ax, CurrentPositionObject)
            % create a group object and group all plots to this "parent" handle
            graphicsHandle = hggroup('Parent', ax);
%             keyVal = ['sagital',this.String,ax];
            keyVal = ax;
%             [base, direcVec] = this.getElectrodePoints;
%             curPosInMM = repmat(CurrentPositionObject.currentImagePosInMmS,1,5);
%             factor = repmat(((curPosInMM - base(1,:))./direcVec(1,:)),3,1);
%             projectionPoint = base + factor.*direcVec;
            curPosInMM = CurrentPositionObject.currentImagePosInMmS;
            base = this.targetPoint;
            direcVec = this.entryPoint - base;
            factor = ((curPosInMM - base(1))./direcVec(1));
            projectionPoint = base + factor.*direcVec;
            
            pointStart = [this.entryPoint(2) this.targetPoint(2)];
            pointEnd = [this.entryPoint(3) this.targetPoint(3)];
            
            projectionPointStart = [this.entryPoint(2) projectionPoint(2)];
            projectionPointEnd = [this.entryPoint(3) projectionPoint(3)];
            
            posLineColor = [0 1 1];
            
            %check wheather we are within entry-target or >taget
            if(projectionPoint(3) < min(this.targetPoint(3),this.entryPoint(3)))%note it starts from left-bottom-corner
                projectionPointStart = [this.targetPoint(2) projectionPoint(2)];
                projectionPointEnd = [this.targetPoint(3) projectionPoint(3)];
                posLineColor = [1 0 0];
            elseif(projectionPoint(3) > max(this.entryPoint(3),this.targetPoint(3)))
                projectionPointStart = [this.entryPoint(2) projectionPoint(2)];
                projectionPointEnd = [this.entryPoint(3) projectionPoint(3)];
                posLineColor = [1 0 0];
            end%if
            
            % Start to Plot what needs to be plotted
            handle.CurrentPositionObject = CurrentPositionObject;
%             this.plotHandleHashMap(1.0) = CurrentPositionObject;
            this.onOrthogonalProjectionPlot(keyVal, graphicsHandle, pointStart, pointEnd, projectionPointStart, projectionPointEnd, posLineColor, handle);
        end%fcn-initPlotSagital
        
        % 'Callback' is called from inherited 2D-plot-methods (see plotableImage) to start
        % plotting
        function onOrthogonalProjectionPlot(this, keyVal, graphicsHandle, pointStart, pointEnd, projectionPointStart, projectionPointEnd, posLineColor, handle)
            if(this.show2DProjection)
                %check if exists a handle for this axes. if yes only change
                %currentPosition line
                if(isKey(this.plotHandleHashMap,double(keyVal))) 
                     %kill old handle (graphic in axes)
                    value = values(this.plotHandleHashMap,{double(keyVal)});
                    try
                        handlesFromHgGroup = get(value{1}.ProjectionHandle, 'Children'); %two-childs: first->curPosLine, second->electrode
                        hgGroup = value{1}.ProjectionHandle; 

                        % update the line
                        handlesFromHgGroup(1).XData = projectionPointStart;
                        handlesFromHgGroup(1).YData = projectionPointEnd;
                    catch
                        
                    end
                    
                else%if not perform whole plot
                    % TODO:add- on these two points for dragging/moving
                    % them and so change the trajectory
                    graphicsHandles.entryPoint = line(pointStart(1),pointEnd(1), 'Marker', ...
                        'o','Color',this.color,'LineStyle','none',...
                        'MarkerSize', this.diameter, 'Parent', graphicsHandle);
                 %   set(graphicsHandles.entryPoint, 'Parent', graphicsHandle);
                    
                    graphicsHandles.targetPoint = line(pointStart(2),pointEnd(2), 'Marker', ...
                        'o','Color',this.color,'LineStyle','none',...
                        'MarkerSize', this.diameter, 'Parent', graphicsHandle);
                  %  set(graphicsHandles.targetPoint, 'Parent', graphicsHandle);
                    
                    graphicsHandles.electrode = line(pointStart,pointEnd, 'Marker', ...
                        this.marker,'Color',this.color,'LineStyle',this.lineStyle,...
                        'MarkerSize', this.diameter, 'Parent', graphicsHandle);
                  %  set(graphicsHandles.electrode, 'Parent', graphicsHandle);
                    
                    graphicsHandles.currentPosLine = line(projectionPointStart,projectionPointEnd, 'Marker', ...
                        this.marker,'Color',posLineColor,'LineStyle','-',...
                        'MarkerSize', this.diameter, 'Parent', graphicsHandle);
                  %  set(graphicsHandles.currentPosLine, 'Parent', graphicsHandle);
                    
                    handle.ProjectionHandle = graphicsHandle;
                    this.plotHandleHashMap(double(keyVal)) = handle;  %speichere hggroup in HashMap
                end%if
                
                set(keyVal, 'DataAspectRatio',[1 1 1], ...
                    'PlotBoxAspectRatioMode','auto');
                set(keyVal, 'Position', [0 0 1 1]); %fill whole axis
                
                % if show2DProjection-Flag ist not set, kill all
           else
%                 this.onHelpDeleteItemsInHashMap({keyVal});
            end%if
        end%fcn-onOrthogonalProjectionPlot
          
        function onUpdatePlots(this, ~, ~)
            this.updatePlot3D;
            this.updatePlot2D;
        end%fcn-onUpdatePlots
        
        %TODO: think about removing this configurable part. not needed
        %anymore
        function panel = getConfigPanel(this, varargin)
%             stepSize = 0.5;
%             if(this.trajectoryLength > 0)
%                 sliderStepSize = [stepSize, stepSize] / this.trajectoryLength;
%             else
%                 error('TestElectrodes:getConfigPanel', 'trajectory lenght have to be greater than zero');
%             end%if
            
            parent = findArg('Parent', varargin{:});
            panel = uiextras.Panel('Parent',parent,'Title','TestElectrode','Padding', 3);
            
            vBox = uiextras.VBox('Parent',double(panel), 'Spacing', 3);
            
            %horizontal box for checkboxes
            checkBoxBox = uiextras.HBox('Parent', vBox, 'Spacing', 3);
            uicontrol('Parent',checkBoxBox,'Style','checkbox', ...
                'Callback', @this.onCheckBoxTestElectrodeClicked, ...
                'String',   '2D-Projection', ...
                'Value', this.show2DProjection);
            
            HBoxBox = uiextras.HBox('Parent', vBox, 'Spacing', 2);
            
            uicontrol('Parent', HBoxBox,'Style','pushbutton',...
                'String', 'transfrom',...
                'Callback', @this.onTransformationClicked);
            
            uicontrol('Parent', HBoxBox,'Style','pushbutton',...
                'String', 'save',...
                'Callback', @this.onSaveClicked);

            
%             uicontrol('Parent', vBox,'Style','slider', ...
%                 'Min', 0, 'Max', this.trajectoryLength, 'SliderStep', [0.1 0.5] ./ this.trajectoryLength, ...
%                 'Value', 0, 'Callback', @this.sliderMoved);
            

%             set(HBoxBox,'sizes',[-1,-1]);
        end%fcn-getConfigPanel
        
        function onCheckBoxTestElectrodeClicked(this, src, ~)
            val = get(src, 'Value');
            if(isequal(val,1))
                this.show2DProjection = 1;
                this.updatePlot2D();
            else
                this.show2DProjection = 0;
                this.onHelpDeleteItemsInHashMap();
            end%if
        end%fcn-onCheckBoxTestElectrodeClicked        
        
        function onTransformationClicked(this, ~ ,~)
            [fileName,pathName] = uigetfile({'*.nii*'},'Select the target System (only RigidMni and RigidIntraCT)');
            [folder, ~, ~] = fileparts(fullfile(pathName, fileName));
            [folderWithoutRefSys, ref] = fileparts(folder);
            if(~(isequal(ref,'rigidMni') || isequal(ref,'rigidIntraCT')))
                HintWindow('At the moment there is only a transformation between RigidMni and RigidIntraCT possible. Please choose a target Image from either one of this Reference-Systems');
                return
            end%if
            [~,patId] = fileparts(folderWithoutRefSys);
                
%             if !=reverse % transfrom from RigidIntraCT -> RigidMni
            reverse = 1 - isequal(ref,'rigidMni');
            points2RigidMni = {this.entryPoint, this.targetPoint};
            transformedPointsInMmfromRigidIntraCtToRigidMni = transformPointsFromRigidIntraCTToRigidMni(patId, points2RigidMni,reverse);
            [entry, target] = transformedPointsInMmfromRigidIntraCtToRigidMni{:};
            this.entryPoint = entry;
            this.targetPoint = target;
        end%fcn-onTransformationClicked
        
        
        function onSaveClicked(this, ~, ~)
           [fileName, pathName] = uiputfile({[this.toString, '.ele']},'Select Path to store the trajectory');
           save(fullfile(pathName,fileName), 'this'); % save as mat-file
        end%fcn-onSaveClicked
        
        function delete(this) %HN:if an testelectrode is deleted, so this method guarantees that all plots will be deleted$added 16.4.14 
            delete@plotableImage(this); % call superclass destructor
            delete@plotable3D(this);    % call superclass destructor
%             delete@Trajectory(this);
        end%   
    end
    methods(Access=public)
        
        function onHelpDeleteItemsInHashMap(this, keyVal)
            if(nargin < 2)
                keyVal = this.plotHandleHashMap.keys();
            end%if
            if(isKey(this.plotHandleHashMap,keyVal)) 
                    value = values(this.plotHandleHashMap,keyVal);
                    for i=1:length(value)
                        try
                            handlesFromHgGroup = get(value{i}.ProjectionHandle, 'Children'); %two-childs: first->curPosLine, second->electrode
                            %                     handlesFromHgGroup = get(value{1}.ProjectionHandle; %two-childs: first->curPosLine, second->electrode
                            delete(handlesFromHgGroup);
                        catch e
                            disp(['TestElectrodes:onHelpDeleteItemsInHashMap:CantGetItem']);
                        end%try/catch
                    end%for
%                     remove(this.plotHandleHashMap, keyVal);
            end%if
        end%fcn-onHelpDeleteItemsInHashMap
        
    end%methods-private
end%class

%%
function test()

elec = TestElectrodes([0; 150; 200], [80;111; 67]); %RAI

elec.intersections = [1 12 2 0 -5; 10 6 6 6 6]; % =[tg1c tg1a tg1p tg1l tg1r; tg2c tg2a tg2p tg2l tg2r]
elec.initPlot3D(gca);

nii = NiftiSegmentationSubvolume('segmentation_swan_fb.nii');
nii.initPlot3D(gca);
nii.pointCloudInMm('SNr+STN_L')


% d = elec.electrodesTargetPoints
% e = elec.electrodesEntryPoints
% tg = elec.getIntersectionPointCloud
% scatter3(d(1,:),d(2,:),d(3,:), 'b*')
% camtarget(d(:,1));
% hold on
% scatter3(e(1,:),e(2,:),e(3,:), 'g*')
% 
% scatter3(tg(1,:),tg(2,:),tg(3,:), 'r')
% parentAxes = gca;
% line([tg(1,1:5); tg(1,6:10)], ...
%      [tg(2,1:5); tg(2,6:10)], ...
%      [tg(3,1:5); tg(3,6:10)], ...
%      'LineWidth', 2.0,...
%      'Color', 'r',...                           
%      'Parent', parentAxes)

end