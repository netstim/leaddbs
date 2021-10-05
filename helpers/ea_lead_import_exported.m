classdef ea_lead_import_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        LeadMigrateUIFigure          matlab.ui.Figure
        RunButton                    matlab.ui.control.Button
        ConverttoDicomusingDropDownLabel  matlab.ui.control.Label
        ConverttoDicomusingDropDown  matlab.ui.control.DropDown
        MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox  matlab.ui.control.CheckBox
        SelectpatientsButton         matlab.ui.control.Button
        ConvertDICOMdatatoBIDScompliantdatasetCheckBox  matlab.ui.control.CheckBox
        SelectoutputdirectoryButton  matlab.ui.control.Button
        BusyLamp                     matlab.ui.control.Lamp
    end

    
    properties (Access = public)
        source_files % Description
        output_files
        dicom_source_file
        input_label
        output_label
        flag
        close_flag
        handles
        options
        subjID
    end
    
    methods (Access = public)
        
        function results = enablerules(app,selpath)
                for selection = 1:length(selpath)
                    [filepath,app.subjID,ext] = fileparts(selpath);
                    app.subjID = erase(app.subjID,'_');
                    files_inside = dir_without_dots(selpath{selection});
                    files_inside = {files_inside.name};
                    %detection function: if the pat has derivatives and
                    % raw data, but not in BIDS format then it should be
                    % migrated.
                    if any(ismember(files_inside,'glanat.nii')) && any(cellfun('isempty',regexp(files_inside,'raw_anat.*.nii'))) && all(cellfun('isempty',regexpi(files_inside,'DICOM')))  && all(cellfun('isempty',regexp(files_inside,'.*.dcm')))
                        %msgbox("Looks like preprocessing has been completed, but no dicom images are present. Conitinuing to migrate into BIDS-compliant dataset..");
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Enable = 'off';
                        app.ConverttoDicomusingDropDown.Enable = 'off';
                        app.ConverttoDicomusingDropDownLabel.Enable = 'off';
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 1;
                        app.flag = 'onlyMigrate';
                        %only raw nifti available, no dicom is available. Migrate
                        %will be done, but only to move & rename the raw files
                    elseif ~any(ismember(files_inside,'glanat.nii')) && any(cellfun('isempty',regexp(files_inside,'raw_anat.*.nii'))) && ~any(ismember(files_inside,'rawdata')) && all(cellfun('isempty',regexpi(files_inside,'DICOM')))
                        %msgbox("Raw Nifti detected without any DICOM images. Continuing to run Migration to generate raw BIDS-compliant dataset only..");
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Enable = 'off';
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value = 0;
                        app.ConverttoDicomusingDropDown.Enable = 'off';
                        app.ConverttoDicomusingDropDownLabel.Enable = 'off';
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 1;
                        app.flag = 'onlyMigrate';
                        %only DICOM files available
                    elseif ~any(endsWith(files_inside,'.nii')) && (any(ismember(files_inside,'DICOM')) || ~any(cellfun('isempty',regexp(files_inside,'.*.dcm'))))
                        %todo: remove support for .dcm files inside
                        %if any(ismember(files_inside,'DICOM')) || ~any(cellfun('isempty',regexp(files_inside,'.*.dcm')))
                        %msgbox("DICOM files detected. Continuing to run dicom to nifti conversion, and nifti to Lead-BIDS conversion.");
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 0;
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Enable = 'off';
                        app.ConverttoDicomusingDropDown.Enable = 'on';
                        app.ConverttoDicomusingDropDownLabel.Enable = 'on';
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value = 1;
                        app.flag = 'MigrateDCMconv';
                        app.dicom_source_file = fullfile(selpath,'DICOM');
                        %BIDS Compliant raw data nifti files are available! in this case,
                        %dcm -> bids conversion should require special
                        %handling. Only the conversion from nifti -> lead bids
                        %will be performed
                    elseif any(ismember(files_inside,'rawdata'))
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 0;
                        app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Enable = 'off';
                        app.ConverttoDicomusingDropDown.Enable = 'on';
                        app.ConverttoDicomusingDropDownLabel.Enable = 'on';
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Enable = 'on';
                        app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value = 1;
                        app.flag = 'onlyDCMconv';
                        %both dicom & derivatives available. Requires special
                        %handling in the migrate code.
                    elseif any(ismember(files_inside,'glanat.nii')) && any(cellfun('isempty',regexp(files_inside,'raw_anat.*.nii')))
                        if any(ismember(files_inside,'DICOM')) || ~any(cellfun('isempty',regexp(files_inside,'.*.dcm')))
                            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 1;
                            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value = 1;
                            app.ConverttoDicomusingDropDown.Enable = 'on';
                            app.ConverttoDicomusingDropDownLabel.Enable = 'on';
                            app.dicom_source_file = fullfile(selpath,'DICOM');
                            app.flag = 'MigrateDCMconv';
                        end
                    else
                        msgbox("Unable to detect the type of processing required, please select options manually..")
                    end
                end
            end
            %selpath now contains a list of input directories.  
        
        function results = op_label(app,selpath)
           if isempty(app.input_label)
                if length(selpath) == 1
                    app.input_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",'Loaded Source Directory(s)!',"Position",[109,189,120,22]);
                    app.SelectpatientsButton.Text = selpath{1};
                else
                   app.input_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",'Loaded Source Directory(s)!',"Position",[109,189,120,22]);
                   app.SelectpatientsButton.Text = ['Loaded multiple patients(' length(selpath), ')']; 
                end
            else
                app.input_label.Label.Text = '';
                if length(selpath) == 1
                    app.input_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",'Loaded Source Directory(s)!',"Position",[37,198,250,22]);
                    app.SelectpatientsButton.Text = selpath{1};
                else
                   app.input_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",'Loaded Source Directory(s)!',"Position",[37,198,250,22]);
                   app.SelectpatientsButton.Text = ['Loaded multiple patients(' num2str(length(selpath)) ')']; 
                end
            end 
        end
        
        function results = reset_rules(app)
            app.ConverttoDicomusingDropDown.Enable = 'on';
            app.ConverttoDicomusingDropDownLabel.Enable = 'on';
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Enable = 'on';
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value = 0;
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Enable = 'on';
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value = 0;
        end
        function [bidsroot,subjID] = generateopparams(app,dest)
            bidsroot = dest;
            subjID = app.subjID;
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, source_files, options, handles)
            app.LeadMigrateUIFigure.Name = 'Lead Import';
            if exist('source_files','var')
               op_label(app,source_files);
               enablerules(app,source_files);
               app.source_files = source_files;
               app.close_flag = 1;
           end
           if exist('handles','var')
               app.handles = handles;
           end
           if exist('options','var')
               app.options = options;
           end
        end

        % Button pushed function: SelectpatientsButton
        function SelectpatientsButtonPushed(app, event)
            app.BusyLamp.Color = 'red';
            %reset everything
            reset_rules(app)
            selpath = ea_uigetdir(pwd);
            op_label(app,selpath);
            enablerules(app,selpath);
            app.source_files = selpath;
            app.BusyLamp.Color = 'green';
        end

        % Value changed function: 
        % ConvertDICOMdatatoBIDScompliantdatasetCheckBox
        function ConvertDICOMdatatoBIDScompliantdatasetCheckBoxValueChanged(app, event)
            app.BusyLamp.Color = 'red';
            value = app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value;
            if value && ~app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'onlyDCMconv';
            elseif ~value && app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'onlyMigrate';
            elseif value && app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'MigrateDCMconv';
            end
            app.BusyLamp.Color = 'green';
        end

        % Value changed function: ConverttoDicomusingDropDown
        function ConverttoDicomusingDropDownValueChanged(app, event)
            app.BusyLamp.Color = 'red';
            value = app.ConverttoDicomusingDropDown.Value;
            
            if app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value && ~app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'onlyDCMconv';
            elseif app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value && app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'MigrateDCMconv';
            elseif ~app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value && app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value
                app.flag = 'onlyMigrate';
            end
            app.BusyLamp.Color = 'green';
        end

        % Button pushed function: RunButton
        function RunButtonPushed(app, event)
            
            
            app.BusyLamp.Color = 'red';
            source = app.source_files;
            if isempty(app.output_files)
               errordlg("No output directory selected, select to continue...",'Output directory error')
            else
                dest = app.output_files;
                
                %for johannes
                convert_using = app.ConverttoDicomusingDropDown.Value;
                if strcmp(app.flag,'onlyMigrate')
                    f = waitbar(0,'Migrating dataset..');
                    waitbar(0.33,f,'Migrating dataset...');
                    if ~isempty(app.dicom_source_file)
                       legacy2bids(source,dest,1,app.dicom_source_file,1);
                    else
                       legacy2bids(source,dest,0,'',0);
                    end
                    waitbar(1,f,'Done!');
                    close(f);
                    % after migration is done, trigger ea_dataset_import.m
                    % without DICOM conversion to sort out .nii.gz files in BIDS
                    % dataset
                    ea_dataset_import(source, [], [], 0)
                    
                end
                %both the migrate and the dcm code here. Migrate runs first to
                %dump the legacy raw dataset into the new folder, "legacy_rawdata"
                if strcmp(app.flag,'MigrateDCMconv')
                    %app.dicom_source: parent dicom,
                    %implement checkpoint if folder is not empty
                    legacy2bids(source,dest,1,app.dicom_source_file,1);
                    ea_dataset_import(source, dest, convert_using, 1);
                end
                % only DICOM conversion, no migration
                if strcmp(app.flag,'onlyDCMconv')
                    ea_dataset_import(source, dest, convert_using, 1);
                end
                app.BusyLamp.Color = 'green';
                ea_busyaction('on',app.handles.leadfigure,'dbs');
                % if running automatically, close as soon as the migration
                % is over.
                if app.close_flag
                   setappdata(app.handles.leadfigure, 'BIDSRoot', dest{1});
                   setappdata(app.handles.leadfigure, 'subjID', app.subjID);
                   close(app.LeadMigrateUIFigure);
                end
            end
        end

        % Callback function
        function MigrateOptionsButtonGroupSelectionChanged(app, event)
           
        end

        % Value changed function: 
        % MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox
        function MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBoxValueChanged(app, event)
            value = app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value;
            if value && ~app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag = 'onlyMigrate';
            elseif value && app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag = 'MigrateDCMconv';
            elseif ~value && app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag =  'onlyDCMconv';
            end
        end

        % Button pushed function: SelectoutputdirectoryButton
        function SelectoutputdirectoryButtonPushed(app, event)
            app.BusyLamp.Color = 'red';
            app.output_files = ea_uigetdir(pwd);
            if length(app.output_files) == 1
                output_message = app.output_files{1};
            else
                output_message = ['Loaded multiple directories(' num2str(length(app.output_files)) ')'];
            end
            if isempty(app.output_label)
                app.output_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",output_message,"Position",[37,290,250,22]);
            else
                app.output_label.Label.Text = '';
                app.output_label.Label = uilabel(app.LeadMigrateUIFigure,"Text",output_message,"Position",[37,290,250,22]);
            end
            
            %op_status.Label = uilabel(app.LeadMigrateUIFigure,"Text",'Loaded Output Directory(s)!',"Position",[420,290,250,22]);
            if app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value && ~app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag = 'onlyMigrate';
            elseif app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value && app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag = 'MigrateDCMconv';
            elseif ~app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Value && app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Value
                app.flag = 'onlyDCMconv';
            end
            if ~isempty(app.source_files) && ~isempty(app.output_files)
                app.RunButton.Enable = 'on';
            end
            app.BusyLamp.Color = 'green';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create LeadMigrateUIFigure and hide until all components are created
            app.LeadMigrateUIFigure = uifigure('Visible', 'off');
            app.LeadMigrateUIFigure.Color = [1 1 1];
            app.LeadMigrateUIFigure.Position = [100 100 624 285];
            app.LeadMigrateUIFigure.Name = 'Lead Migrate';

            % Create RunButton
            app.RunButton = uibutton(app.LeadMigrateUIFigure, 'push');
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.BackgroundColor = [1 1 1];
            app.RunButton.Enable = 'off';
            app.RunButton.Position = [513 42 100 22];
            app.RunButton.Text = 'Run';

            % Create ConverttoDicomusingDropDownLabel
            app.ConverttoDicomusingDropDownLabel = uilabel(app.LeadMigrateUIFigure);
            app.ConverttoDicomusingDropDownLabel.HorizontalAlignment = 'right';
            app.ConverttoDicomusingDropDownLabel.Position = [331 102 139 22];
            app.ConverttoDicomusingDropDownLabel.Text = 'Convert to Dicom using..';

            % Create ConverttoDicomusingDropDown
            app.ConverttoDicomusingDropDown = uidropdown(app.LeadMigrateUIFigure);
            app.ConverttoDicomusingDropDown.Items = {'dcm2niix', 'dicm2nii', 'SPM'};
            app.ConverttoDicomusingDropDown.ValueChangedFcn = createCallbackFcn(app, @ConverttoDicomusingDropDownValueChanged, true);
            app.ConverttoDicomusingDropDown.BackgroundColor = [1 1 1];
            app.ConverttoDicomusingDropDown.Position = [487 102 126 22];
            app.ConverttoDicomusingDropDown.Value = 'dcm2niix';

            % Create MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox = uicheckbox(app.LeadMigrateUIFigure);
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.ValueChangedFcn = createCallbackFcn(app, @MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBoxValueChanged, true);
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Text = 'Migrate legacy Lead-DBS dataset to BIDS dataset';
            app.MigratelegacyLeadDBSdatasettoBIDSdatasetCheckBox.Position = [29 154 296 22];

            % Create SelectpatientsButton
            app.SelectpatientsButton = uibutton(app.LeadMigrateUIFigure, 'push');
            app.SelectpatientsButton.ButtonPushedFcn = createCallbackFcn(app, @SelectpatientsButtonPushed, true);
            app.SelectpatientsButton.BackgroundColor = [1 1 1];
            app.SelectpatientsButton.Position = [37 219 250 22];
            app.SelectpatientsButton.Text = 'Select patient(s)';

            % Create ConvertDICOMdatatoBIDScompliantdatasetCheckBox
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox = uicheckbox(app.LeadMigrateUIFigure);
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.ValueChangedFcn = createCallbackFcn(app, @ConvertDICOMdatatoBIDScompliantdatasetCheckBoxValueChanged, true);
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Text = 'Convert DICOM data to BIDS compliant dataset';
            app.ConvertDICOMdatatoBIDScompliantdatasetCheckBox.Position = [28 102 297 22];

            % Create SelectoutputdirectoryButton
            app.SelectoutputdirectoryButton = uibutton(app.LeadMigrateUIFigure, 'push');
            app.SelectoutputdirectoryButton.ButtonPushedFcn = createCallbackFcn(app, @SelectoutputdirectoryButtonPushed, true);
            app.SelectoutputdirectoryButton.BackgroundColor = [1 1 1];
            app.SelectoutputdirectoryButton.Position = [363 219 250 22];
            app.SelectoutputdirectoryButton.Text = 'Select output directory';

            % Create BusyLamp
            app.BusyLamp = uilamp(app.LeadMigrateUIFigure);
            app.BusyLamp.Position = [593 255 20 20];

            % Show the figure after all components are created
            app.LeadMigrateUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = ea_lead_import_exported(varargin)

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.LeadMigrateUIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.LeadMigrateUIFigure)
        end
    end
end