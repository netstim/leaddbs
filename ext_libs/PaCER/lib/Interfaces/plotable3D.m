%% plotable3D - Abstract class to be inherited by classes providing 3D plot functions
% plotable3D matches automatic plot redrawing on changes of desired
% properties (by adding approriate listeners there)
%
% Andreas Husch, Florian Bernard
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, fbernardpi@gmail.com

classdef plotable3D < handle
    properties (Access = public)
        plotHandleHashMap3D = [];       
    end
    methods(Abstract)
        handles = initPlot3D(axes)
        string = toString(axes);
    end
    methods
        function this = plotable3D()
            this.plotHandleHashMap3D = containers.Map('KeyType', 'double', 'ValueType', 'double');
        end
        function updatePlot3D(this, ~, ~) % child class should assign listerns to appropiate properties to call this update function
            handles = this.plotHandleHashMap3D.values;
            keys = this.plotHandleHashMap3D.keys; %==axes
            for i=1:this.plotHandleHashMap3D.Count
                this.plotHandleHashMap3D.remove(keys{i});
                delete(handles{i}); %kill the graphic objects

                if(ishandle(keys{i})) % check that the axes still exists
                    this.initPlot3D(keys{i}); %replot
                end
            end
        end
        function delete(this) % destructor
            handles = this.plotHandleHashMap3D.values;
            for i=1:this.plotHandleHashMap3D.Count
                if(ishandle(handles{i}))
                    delete(handles{i}); %kill the graphic objects
                end
            end
            keys = this.plotHandleHashMap3D.keys;
            remove(this.plotHandleHashMap3D, keys);
        end
        function handles = plot3D(this)
            handles = this.initPlot3D(gca);
        end
    end
end