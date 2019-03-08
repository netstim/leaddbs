%% plotableImage - Abstract Class to inherit by objects providing 2D plot methods
%
% Andreas Husch, Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, fbernardpi@gmail.com

classdef (Abstract) plotableImage < handle
    properties (Access = public)
        plotHandleHashMap = [];
    end
    methods(Abstract)
        initPlotAxial(axes, CurrentPositionObject)
        initPlotCoronal(axes, CurrentPositionObject)
        initPlotSagital(axes, CurrentPositionObject)
        toString(axes);
    end
    methods % abstract class including concretem methods!
        function this = plotableImage()
            this.plotHandleHashMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
                  
        function delete(this) % destructor
            handles = this.plotHandleHashMap.values;
            for i=1:this.plotHandleHashMap.Count
                if(ishandle(handles{i}.ProjectionHandle))
                    delete(handles{i}.ProjectionHandle); %kill the graphic objects
                end
            end
            keys = this.plotHandleHashMap.keys;
            remove(this.plotHandleHashMap, keys);
        end%fcn-delete
        
        function updatePlot2D(this, ~, ~) % child class should assign listerns to appropiate properties to call this update function
            handles = this.plotHandleHashMap.values;
            keys = this.plotHandleHashMap.keys;
            
            for i=1:this.plotHandleHashMap.Count 
                CurrentPositionObject = handles{i}.CurrentPositionObject;
                if(ishandle(handles{i}.ProjectionHandle))
                    delete(handles{i}.ProjectionHandle); %kill the graphic objects
                end
                this.plotHandleHashMap.remove(keys{i});
                
                if(ishandle(keys{i})) % check that the axes still exists
                    axUserData = get(keys{i},'UserData');
                    if(isempty(axUserData)) % HA: Workaround bug in oblique View!!
                        disp('plotableImage: UserData of Axes not set');
                        return;
                    end
                    switch axUserData.VisuType
                        case 'axial'
                            this.initPlotAxial(keys{i}, CurrentPositionObject);
                        case 'coronal'
                            this.initPlotCoronal(keys{i}, CurrentPositionObject);
                        case 'sagital'
                            this.initPlotSagital(keys{i}, CurrentPositionObject);
                        case 'ObliqueView'
                            disp('plotableImage:updatePlot2D: missing case');
                        otherwise %do nothing
                    end
                end
                
            end
        end
        
        % should be called if 2d-rpojections should be killed
        function onHelpDeleteProjections(this, keyVal)
            if(nargin < 2)
                keyVal = this.plotHandleHashMap.keys();
            else
                keyVal = {double(keyVal)};
            end%if
            if(isKey(this.plotHandleHashMap,keyVal)) 
                    value = this.plotHandleHashMap.values(keyVal);
                    for i=1:length(value)
                        try
%                             if(ishandle(value{i}.ProjectionHandle) || ishghandle(value{i}.ProjectionHandle)) %two-childs: first->curPosLine, second->electrode
                            %                     handlesFromHgGroup = get(value{1}.ProjectionHandle; %two-childs: first->curPosLine, second->electrode
                            delete(value{i}.ProjectionHandle);
                        catch
                            disp('plotbaleImage:onHelpDeleteProjections:CantDeleteItem');
                        end
                    end
                    this.plotHandleHashMap.remove(keyVal);
            end
        end
    end
end%
