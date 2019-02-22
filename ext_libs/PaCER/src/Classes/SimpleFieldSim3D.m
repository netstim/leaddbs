%% Visualisation of the Mädler/Coenen VTA Approximation
% see Mädler, B., and V. A. Coenen. "Explaining clinical effects of deep brain stimulation 
% through simplified target-specific modeling of the volume of activated tissue." 
% American Journal of Neuroradiology 33.6 (2012): 1072-1080.
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2013 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
classdef SimpleFieldSim3D < plotable3D & configurable & handle
    properties (SetAccess = public, GetAccess = public)
        trajectoryObject = [];
        color = [1 0 0];    %default red
        string = 'SimpleFieldSim';
        graphicsHandles = [];
    end
    properties (SetAccess = public, GetAccess = public, SetObservable=true)
        impedance = 800; %[Ohm]
        voltage = 1.0; %[V]
    end
    properties (Access = protected)
        sphere_x = []; % a normalized sphere used as base
        sphere_y = [];
        sphere_z = [];
    end
    
    methods
        function this = SimpleFieldSim3D(trajectoryObject)
            this.trajectoryObject = trajectoryObject;
            % init normalized sphere
            [this.sphere_x,this.sphere_y,this.sphere_z] = sphere;
            addlistener(this, 'voltage', 'PostSet', @this.updatePlot3D);
            addlistener(this, 'impedance', 'PostSet', @this.updatePlot3D);
            if(isa(this.trajectoryObject, 'MetaTrajectory'))
                addlistener(this.trajectoryObject, 'trajectoryChanged', @this.updatePlot3D); %inherited method
            end
        end
        
        function str = toString(this)
            str = this.string;
        end
        
        function graphicsHandle = initPlot3D(this, ax)
            set(ax, 'NextPlot', 'add');
            if(isa(this.trajectoryObject, 'TestElectrodes'))
                targetpoint = this.trajectoryObject.electrodesTargetPoints(:,this.selectedElectrodeContact); %NOT contact but the selected trajectory (out of the five)
            elseif(isa(this.trajectoryObject, 'PolynomialElectrodeModel'))
                targetpoint = this.trajectoryObject.activeContactPoint; 
            else
                targetpoint = this.trajectoryObject.targetPoint;
            end
            r = this.calcFieldradius(this.impedance, this.voltage);
            x = this.sphere_x .* r + targetpoint(1);
            y = this.sphere_y .* r + targetpoint(2);
            z = this.sphere_z .* r + targetpoint(3);

            graphicsHandle = surf(x,y,z, 'EdgeColor', 'none', ...
                'FaceAlpha', 0.6, 'FaceColor',this.color,...
                'FaceLighting', 'gouraud', 'Parent', ax);
            this.plotHandleHashMap3D(double(ax)) = double(graphicsHandle);
            
            daspect([1 1 1]);
            % add light if none is present
            if(isempty(findobj(ax, 'Type', 'light')))
                camlight headlight;
            end
        end

        function panel = getConfigPanel(this, varargin)
            panel = SimpleFieldSim3DViewElement(this, varargin{:});
        end
     
        function fig = getConfigWindow(this)
            scrsz = get(groot,'ScreenSize');
            fig = figure('Name', 'SimpleFieldSim3D Config', 'Position',[scrsz(3)/2 scrsz(4)/2 500 150], 'MenuBar', 'none');
            this.getConfigPanel();
        end
    end   
    
    methods (Static)
        function r = calcFieldradius(impedance, voltage)
            k1 = -1.0473;
            k3 = 0.2786;
            k4 = 0.0009856;
            
            % calculate isoline distance [mm] from center of electrode. note the electrode type!
            fr = @(impedance, voltage)(-(k4 * impedance -sqrt(k4^2*impedance^2 + 2*k1*k4*impedance + k1^2 + 4*k3*voltage) + k1) / (2 * k3));
            r = fr(impedance, voltage);
        end
    end
end