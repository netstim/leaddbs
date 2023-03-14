function [hObject, eventdata, handles] = ea_convertGUIDEArguments(app, event)
    %CONVERTTOGUIDECALLBACKARGUMENTS Converts App Designer callback
    %arguments to GUIDE-style callback arguments.
    %   [HOBJECT, EVENTDATA, HANDLES] = CONVERTTOGUIDECALLBACKARGUMENTS(APP, EVENT)
    %   is automatically added to the beginning of each migrated app
    %   callback by the GUIDE to App Designer Migration Tool to enable
    %   GUIDE-style callback code to be functional in an App Designer
    %   callback. It creates GUIDE-style callback arguments HOBJECT,
    %   EVENTDATA, and HANDLES from App Designer callback arguments APP and
    %   EVENT.
    %
    %   HOBJECT is the handle of the object whose callback is executing.
    %   For components that were of type uicontrol before being migrated to
    %   App Designer, HOBJECT will be a uicontrol-like object that has the
    %   same properties as a uicontrol to enable most uicontrol-style code
    %   to function in an App Designer callback.
    %
    %   EVENTDATA is usually empty but can be a structure containing
    %   specific information about the callback event.
    %
    %   HANDLES is a structure containing all the child components of the
    %   app's figure that have a value defined for their Tag property.
    %   Child components of the figure that were of type uicontrol before
    %   being migrated to App Designer, will be a uicontrol-like object
    %   that has the same properties as a uicontrol to enable most
    %   uicontrol-style code to function in an App Designer callback.

    % Copyright 2019-2020 The MathWorks, Inc.
    narginchk(1, 2);

    appFigure = appdesigner.internal.service.AppManagementService.getFigure(app);

    if isempty(appFigure)
        % Some callbacks (i.e. SizeChangedFcn) are called before App
        % Designer has associated the app and its figure.  In that case,
        % look at the figure ancestor of the event source.
        appFigure = ancestor(event.Source, 'figure');
    end

    handles = guidata(appFigure);

    if isempty(handles)
        handles = localCreateHandles(appFigure);
        guidata(appFigure, handles);
    end
    handles.prod='dbs';
    handles.callingfunction='lead_dbs';
    handles.leadfigure=app.leadfigure;
    % StartupFcn case - there is no event or eventdata, and the hObject
    % should be the app's figure.
    if nargin == 1
        hObject = appFigure;
        eventdata = [];
        return;
    end

    % Retrieve the source component for the event.  In most cases, this
    % will be the value for hObject.
    sourceComponent = event.Source;

    % By default, the eventdata is the same as the event.
    eventdata = event;

    % For buttongroup SelectionChangedEvent, hObject in GUIDE was set to
    % the UIControl RadioButton that is newly selected.  This newly
    % selected radio button is stored in the NewValue property of the event
    %
    % Use the EventName property to allow testing of this function
    if strcmp(event.EventName, 'SelectionChanged') && strcmp(sourceComponent.Type, 'uibuttongroup')
        sourceComponent = event.NewValue;

        if isprop(event.NewValue, 'CodeAdapter')
            newValue = event.NewValue.CodeAdapter;
        else
            newValue = event.NewValue;
        end

        if isprop(event.OldValue, 'CodeAdapter')
            oldValue = event.OldValue.CodeAdapter;
        else
            oldValue = event.OldValue;
        end

        % Create a mock event with code adapters for the buttons, if needed
        eventdata = struct('EventName', 'SelectionChanged', 'NewValue', newValue, 'OldValue', oldValue, 'Source', event.Source);
    end

    % Retrieve the CodeAdapter for this object if it exists, otherwise just
    % return the standard App Designer component.
    if isprop(sourceComponent, 'CodeAdapter')
        hObject = sourceComponent.CodeAdapter;
    else
        hObject = sourceComponent;
    end
end

function handles = localCreateHandles(appFigure)
    % Create a handles struct with all components with their Tag property
    % set.

    % Find all components with Tag not empty.
    migratedComponents = findall(appFigure, '-not','Tag','');
    handles = struct();

    % For all components, assign either the component or its code adapter
    % to the handles struct.  The struct field name is the Tag property.
    for idx = 1:length(migratedComponents)
        component = migratedComponents(idx);
        tag = component.Tag;

        if isprop(component, 'CodeAdapter')
            component = component.CodeAdapter;
        elseif isprop(component, 'Type') 
            switch component.Type
                case {'uibutton', 'uistatebutton', 'uiradiobutton', 'uieditfield', 'uitextarea', 'uislider', 'uilistbox', 'uidropdown', 'uicheckbox', 'uilabel'}
                    component = appdesigner.appmigration.UIControlPropertiesConverter(component);
                case 'uibuttongroup'
                    component = appdesigner.appmigration.ButtonGroupPropertiesConverter(component);
                otherwise
                    % do nothing, the component does not need an adapter
            end
        end

        handles.(tag) = component;
    end
end
