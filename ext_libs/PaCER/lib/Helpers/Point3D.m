%% Point3D - object representing a Point in 3D spaec
%  "no frills" Version 
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicine
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
classdef Point3D < plotable3D  & handle
  
    properties (SetAccess = public, GetAccess = public) 
        point = NaN(3,1);   %x,y,z coordinate in TODO mm
        marker = 'o';       %default circle
        color = [0 1 0];    %default green
        diameter = 3;       % groesse
        string = 'Point3D';
    end
    methods       
        function this = Point3D(varargin)
            if(nargin==1)
                this.point = reshape(varargin{1},3,1);
            end
            
        end
        
        function str = toString(this)
            str = this.string;
        end
        
        function graphicsHandle = initPlot3D(this, ax) %HA: sollte man hier nicht generell ein handle (bzw. eine handle menge) zurückliefern!?
            graphicsHandle= line(this.point(1),this.point(2),this.point(3), 'Marker', ...
                this.marker,'Color',this.color,'LineStyle','none',...
                'MarkerSize', this.diameter, 'Parent', ax);
            this.plotHandleHashMap3D(ax) = graphicsHandle;
        end
        
        function delete(this)
            delete@plotable3D(this);
        end   
    end
    
    methods(Access=private)      
        function onHelpToDeleteItemsFromHashMap(this, keyVal)
            value = values(this.plotHandleHashMap,keyVal);
            delete(value{1}.ProjectionHandle);
            remove(this.plotHandleHashMap, keyVal);
        end
    end
end