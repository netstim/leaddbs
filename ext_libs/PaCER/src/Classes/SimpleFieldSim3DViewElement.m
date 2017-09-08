%% SimpleFieldSim3DViewElement - A panel that represents a View Element
%  to configure a SimpleFieldSim3D Graphics Object
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2013 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
classdef SimpleFieldSim3DViewElement < uiextras.Panel
    properties (Access = protected)
        sldImpedance = [];
        edtImpedance = [];
        sldVoltage = [];
        edtVoltage = [];
        rdoSelectElectrodeContact = [];
        SimpleFieldSim3DObject = SimpleFieldSim3D.empty();
    end
    methods
        function this = SimpleFieldSim3DViewElement(simpleFieldSim3DObject, varargin)
            if(~isa(simpleFieldSim3DObject, 'SimpleFieldSim3D'))
                error('SimpleFieldSim3DViewElement:constructor:invalidParameter', 'First parameter has to be a SimpleFieldSim3D Object');
            end
            this = this@uiextras.Panel(varargin{:});
            this.SimpleFieldSim3DObject = simpleFieldSim3DObject;
            
            %% GUI Elements
            set(this,'TitlePosition','centertop','Padding', 5);
            mainVBox = uiextras.VBox('Parent', this);
            
            %%%%% Sliders for Impedance and Voltage %%%%%
            hBox = uiextras.HBox('Parent',mainVBox);
            uicontrol('Style', 'text', 'Parent', uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'left', 'String', 'Impedance:');
            %min
            uicontrol('Style', 'text', 'Parent', uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'right', 'String', '100');
            this.sldImpedance = uicontrol('Style', 'slider', 'Parent',uiextras.HButtonBox('Parent',hBox, 'ButtonSize', [250 18]),...
                'Min', 100, 'Max', 1400, 'SliderStep', [1/1300 10/1300], 'Value',  this.SimpleFieldSim3DObject.impedance, 'BackgroundColor', 'w',...
                'Callback', @this.onSliderChanged);
            %max
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'left', 'String', '1400');
            this.edtImpedance = uicontrol('Style', 'edit', 'Parent',uiextras.HButtonBox('Parent',hBox),...
                'String', num2str(this.SimpleFieldSim3DObject.impedance), 'BackgroundColor', 'w',...
                'HorizontalAlignment', 'left',...
                'Callback', @this.onEdtChanged);
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'left', 'String', '[Ohm]');
            set(hBox, 'Widths', [80 30 -1 30 40 40])

            hBox = uiextras.HBox('Parent',mainVBox);
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), 'String', 'Voltage:')
            %min
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'right', 'String', '0.5');
            this.sldVoltage = uicontrol('Style', 'slider', 'Parent',uiextras.HButtonBox('Parent',hBox, 'ButtonSize', [250 18]),...
                'Min', 0.5, 'Max', 10,  'SliderStep', [0.1/9.5 0.5/9.5], 'Value',  this.SimpleFieldSim3DObject.voltage, 'BackgroundColor', 'w',...
                'Callback', @this.onSliderChanged);
            %max
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'left', 'String', '10');
            this.edtVoltage = uicontrol('Style', 'edit', 'Parent',uiextras.HButtonBox('Parent',hBox),...
                'String', num2str(this.SimpleFieldSim3DObject.voltage), 'BackgroundColor', 'w',...
                'HorizontalAlignment', 'left',...
                'Callback', @this.onEdtChanged);
            uicontrol('Style', 'text', 'Parent',uiextras.HButtonBox('Parent',hBox), ...
                'HorizontalAlignment', 'left', 'String', '[V]');
            set(hBox, 'Widths', [80 30 -1 30 40 40])
        end
    end
    methods (Access = protected)
        function onSliderChanged(this, ~, ~)
            this.SimpleFieldSim3DObject.voltage = (get(this.sldVoltage, 'Value'));
            this.SimpleFieldSim3DObject.impedance = (get(this.sldImpedance, 'Value'));
            set(this.edtVoltage, 'String', num2str(this.SimpleFieldSim3DObject.voltage));
            set(this.edtImpedance, 'String', num2str(this.SimpleFieldSim3DObject.impedance));
        end
        function onEdtChanged(this, ~, ~)
            this.SimpleFieldSim3DObject.voltage = str2num(get(this.edtVoltage, 'String')); %#ok<ST2NM>
            this.SimpleFieldSim3DObject.impedance = str2num(get(this.edtImpedance, 'String')); %#ok<ST2NM>
            set(this.sldVoltage, 'Value', this.SimpleFieldSim3DObject.voltage);
            set(this.sldImpedance, 'Value', this.SimpleFieldSim3DObject.impedance);
        end
        function onSelectedElectrodeContactChanged(this, ~, ~)
            this.SimpleFieldSim3DObject.selectedElectrodeContact = ...
                this.rdoSelectElectrodeContact.SelectedChild;
        end
    end
end