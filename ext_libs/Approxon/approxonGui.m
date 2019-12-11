function h = approxonGui()

gui = createInterface();
gui.model = load('activation_model_3v.mat'); % loads PW, D, T for 3V 

% Now update the GUI with the current data
redrawDemo();
h = gui.Window;
movegui(h, 'center')
set(h, 'Visible', 'on')


%-------------------------------------------------------------------------%
    function gui = createInterface()
        % Create the user interface for the application and return a
        % structure of handles for global use.
        gui = struct();
        gui.pointHandle = [];
        gui.PW = [];
        gui.D = [];
        % Open a window and add some menus
        gui.Window = figure( ...
            'Name', 'Approximate T by PW, D', ...
            'NumberTitle', 'off', ...
            'MenuBar', 'none', ...
            'Toolbar', 'none', ...
            'HandleVisibility', 'off', ...
            'Position', [100, 100, 600, 500], ...
            'Visible', 'off');
        
        % Arrange the main interface
        mainLayout = uix.VBoxFlex( 'Parent', gui.Window, 'Spacing', 3 );
        
        % + Create the panels
        gui.ViewPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Model Plot          ', ...
            'HelpFcn', @onDemoHelp );
        gui.ViewContainer = uicontainer( ...
            'Parent', gui.ViewPanel );      
        controlPanel = uix.BoxPanel( ...
            'Parent', mainLayout, ...
            'Title', 'Specify PW and D           ' );

        % + Adjust the main layout
        set( mainLayout, 'Heights', [-1,80]  );
        
        controlVBox = uix.VBox( 'Parent', controlPanel, ...
            'Padding', 3, 'Spacing', 3);
        
        % + Create the controls
        controlLayout = uix.HButtonBox( 'Parent', controlVBox, ...
            'Padding', 3, 'Spacing', 3, 'ButtonSize', [100 20], 'VerticalAlignment', 'middle');
        
        gui.TxtD = uicontrol( 'Style', 'text', ...
            'Parent', controlLayout, ...
            'String', 'Axon Diameter [um]:');
        
        gui.EditD = uicontrol( 'Style', 'edit', ...
            'BackgroundColor', 'w', ...
            'Parent', controlLayout, ...
            'String', '3.5', ...
            'Value', 3.5, ...
            'Callback', @onListSelection);
        
        gui.TxtPW = uicontrol( 'Style', 'text', ...
            'Parent', controlLayout, ...
            'String', 'Pulse Width [us]:');
        
        
        gui.EditPW = uicontrol( 'Style', 'edit', ...
            'BackgroundColor', 'w', ...
            'Parent', controlLayout, ...
            'String', '60', ...
            'Value', 60, ...
            'Callback', @onListSelection);
        
        uicontrol( 'Style', 'text', ...
            'Parent', controlLayout, ...
            'String', 'Est. Threshold:');
                
        gui.TxtT = uicontrol( 'Style', 'edit', ...
            'BackgroundColor', [0.5 0.5 0.5], ...
            'Parent', controlLayout, ...
            'String', '60', ...
            'Value', 60, ...
            'Callback', @onListSelection);
        
        
        % ----
        submitBtnBox = uix.HButtonBox( 'Parent', controlVBox, ...
            'Padding', 3, 'Spacing', 3, 'ButtonSize', [100 20], 'VerticalAlignment', 'middle');
        
        gui.SubmitBtn = uicontrol( 'Style', 'pushbutton', ...
            'Parent', submitBtnBox, ...
            'String', 'Save and Return', ...
            'Callback', @onExit);
        
        set( controlVBox, 'Height', [-1 -1] ); % Make the list fill the space
        
        % + Create the view
        p = gui.ViewContainer;
        gui.ViewAxes = axes( 'Parent', p );
        
        
    end % createInterface

%-------------------------------------------------------------------------%
    function redrawDemo()
        % Draw a demo into the axes provided
        
        if isempty ( gui.pointHandle )
                delete( gui.ViewAxes );
                gui.ViewAxes = axes('Parent', gui.ViewContainer);
                gui.pointHandle = plotModelAndGeneralHeuristic(gui.ViewAxes); %#ok<NASGU>
        end
        gui.PW = str2num(gui.EditPW.String); %#ok<ST2NM>
        gui.D = str2num(gui.EditD.String); %#ok<ST2NM>
        assignin('base','PW', gui.PW);   % thats an ugly hack to get our values back, but anyway.. 
        assignin('base', 'D', gui.D);
        
        gui.thresh = str2double(sprintf('%3.2f',  gui.model.activation_model_3v(gui.PW,gui.D)));
        gui.TxtT.String = num2str( gui.thresh);
        delete(gui.pointHandle);
        setappdata(gui.Window, 'thresh', gui.thresh);

        gui.pointHandle = plot3(gui.PW,gui.D,gui.thresh, 'o', 'MarkerSize', 20, ....
         'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r', 'Parent', gui.ViewAxes, 'DisplayName', 'Choosen Parameters'); %#ok<NASGU> %

%         % Some demos create their own figure. Others don't.
%         fcnName = data.DemoFunctions{data.SelectedDemo};
%         switch upper( fcnName )
%             case 'LOGO'
%                 % These demos open their own windows
%                 evalin( 'base', fcnName );
%                 gui.ViewAxes = gca();
%                 fig = gcf();
%                 set( fig, 'Visible', 'off' );
%                 
%             otherwise
%                 % These demos need a window opening
%                 fig = figure( 'Visible', 'off' );
%                 evalin( 'base', fcnName );
%                 gui.ViewAxes = gca();
%         end
%         % Now copy the axes from the demo into our window and restore its
%         % state.
          cmap = colormap( gui.ViewAxes );
          set( gui.ViewAxes, 'Parent', gui.ViewContainer );
          colormap( gui.ViewAxes, cmap );
%         rotate3d( gui.ViewAxes, 'on' );
%         % Get rid of the demo figure
%         close( fig );
    end % redrawDemo

%-------------------------------------------------------------------------%
    function onListSelection( ~, ~ )
        % User selected a demo from the list - update "data" and refresh
        %data.SelectedDemo = get( src, 'Value' );
        redrawDemo();
    end % onListSelection

%-------------------------------------------------------------------------%
    function onExit( ~, ~ )
        % User wants to quit out of the application
        uiresume(gui.Window);
    end % onExit

end % EOF