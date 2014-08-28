function [hPropsPane,parameters] = ea_propertiesGUI(hParent, parameters)
% propertiesGUI displays formatted editable list of properties
%
% Syntax:
%    [hPropsPane,parameters] = propertiesGUI(hParent, parameters)
%
% Description:
%    propertiesGUI processes a list of data properties and displays
%    them in a GUI table, where each parameter value has a unique
%    associated editor.
%
%    propertiesGUI by itself, with no input parameters, displays a demo
%
%    By default, propertiesGUI identifies and processes the following
%    field types: signed, unsigned, float, file, folder, text or string,
%    color, IPAddress, password, date, boolean, cell-array, numeric array,
%    font, struct and class object.
%
% Inputs:
%    hParent - optional handle of a parent GUI container (figure/uipanel
%              /uitab) in which the properties table will appear.
%              If missing or empty or 0, the table will be shown in a
%              new modal dialog window; otherwise it will be embedded
%              in the parent container.
%
%    parameters - struct or object with data fields. The fields are
%              processed separately to determine their corresponding cell
%              editor. If parameters is not specified, then the global
%              test_data will be used. If test_data is also empty, then
%              a demo of several different data types will be used.
%
% Outputs:
%    hPropsPane - handle of the properties panel widget, which can be
%              customized to display field descriptions, toolbar, etc.
%
%    parameters - the resulting (possibly-updated) parameters struct.
%              Naturally, this is only relevant in case of a modal dialog.
%
%    (global test_data) - this global variable is updated internally when
%              the <OK> button is clicked. It is meant to enable easy data
%              passing between the properties GUI and other application
%              component. Using global vars is generally discouraged as
%              bad programming, but it simplifies component interaction.
%
% Customization:
%    This utility is meant to be used either as stand-alone, or as a
%    template for customization. For example, you can attach a unique
%    description to each property that will be shown in an internal
%    sub-panel: see the customizePropertyPane() and preparePropsList()
%    sub-functions.
%
%    When passing the properties in an input parameters struct, the
%    utility automatically inspects each struct field and assigns a
%    corresponding cell-editor with no description and a field label
%    that reflects the field name. The properties are automatically
%    set as modifiable (editable) and assigned a default callback
%    function (propUpdatedCallback() sub-function).
%    See the demoParameters() sub-function for some examples.
%
%    You can have specific control over each property's description,
%    label,  editability, cell-editor and callback function. See the
%    preparePropsList() sub-functions for some examples. You can add
%    additional cell-editors/renderers in the newProperty() sub-function.
%
%    You can place specific control over the acceptable property values
%    by entering custom code into the checkProp() sub-function.
%
% Future development:
%    1. Improve the editor for time, numeric and cell arrays
%    2. Enable more control over appearance and functionality via 
%       propertiesGUI's input parameters
%    3. Add additional built-in cell editors/renderers: slider, point,
%       rectangle (=position), ...
%
% Example:
%    propertiesGUI;   % displays the demo
%
%    params.name   = 'Yair';
%    params.age    = uint8(41);
%    params.folder = pwd;
%    params.date   = now;
%    params.size.width  = 10;
%    params.size.height = 20;
%    [hPropsPane, params] = propertiesGUI(params);
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% Warning:
%    This code heavily relies on undocumented and unsupported Matlab
%    functionality. It works on Matlab 7+, but use at your own risk!
%
%    A technical description of the implementation can be found at:
%    http://undocumentedmatlab.com/blog/propertiesGUI/
%    http://undocumentedmatlab.com/blog/jide-property-grids/
%    http://undocumentedmatlab.com/blog/advanced-jide-property-grids/
%
% Change log:
%    2013-12-24: Fixes for R2013b & R2014a; added support for Font property
%    2013-04-23: Handled multi-dimensional arrays
%    2013-04-23: Fixed case of empty ([]) data, handled class objects & numeric/cell arrays, fixed error reported by Andrew Ness
%    2013-01-26: Updated help section
%    2012-11-07: Minor fix for file/folder properties
%    2012-11-07: Accept any object having properties/fields as input parameter; support multi-level properties
%    2012-10-31: First version posted on <a href="http://www.mathworks.com/matlabcentral/fileexchange/authors/27420">MathWorks File Exchange</a>
%
% See also:
%    inspect, uiinspect (#17935 on the MathWorks File Exchange)

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.09 $  $Date: 2013/12/24 14:33:21 $

  % Get the initial data
  global test_data
  if nargin < 2
      try
          isObj = nargin==1;
          [hasProps,isHG] = hasProperties(hParent);
          isObj = isObj && hasProps && ~isHG;
      catch
          % ignore - maybe nargin==0, so no hParent is available
      end
      if isObj
          parameters = hParent;
          hParent = [];
      else
          parameters = test_data;  % comment this if you do not want persistent parameters
          if isempty(parameters)
              % demo mode
              parameters = demoParameters;
          end
      end
  end

  % Accept any object having data fields/properties
  try
      parameters = get(parameters);
  catch
      oldWarn = warning('off','MATLAB:structOnObject');
      parameters = struct(parameters);
      warning(oldWarn);
  end
  
  % Init JIDE
  com.mathworks.mwswing.MJUtilities.initJIDE;

  % Prepare the list of properties
  oldWarn = warning('off','MATLAB:hg:JavaSetHGProperty');
  warning off MATLAB:hg:PossibleDeprecatedJavaSetHGProperty
  isEditable = true; %=nargin < 1;
  propsList = preparePropsList(parameters, isEditable);
  
  % Create a mapping propName => prop
  propsHash = java.util.Hashtable;
  propsArray = propsList.toArray();
  for propsIdx = 1 : length(propsArray)
      thisProp = propsArray(propsIdx);
      propName = getPropName(thisProp);
      propsHash.put(propName, thisProp);
  end
  warning(oldWarn);

  % Prepare a properties table that contains the list of properties
  model = javaObjectEDT(com.jidesoft.grid.PropertyTableModel(propsList));
  model.expandAll();

  % Prepare the properties table (grid)
  grid = javaObjectEDT(com.jidesoft.grid.PropertyTable(model));
  grid.setShowNonEditable(grid.SHOW_NONEDITABLE_BOTH_NAME_VALUE);
  %set(handle(grid.getSelectionModel,'CallbackProperties'), 'ValueChangedCallback', @propSelectedCallback);
  com.jidesoft.grid.TableUtils.autoResizeAllColumns(grid);
  %com.jidesoft.grid.TableUtils.autoResizeAllRows(grid);
  grid.setRowHeight(19);  % default=16; autoResizeAllRows=20 - we need something in between

  % Auto-end editing upon focus loss
  grid.putClientProperty('terminateEditOnFocusLost',true);

  % If no parent (or the root) was specified
  if nargin < 1 || isempty(hParent) || isequal(hParent,0)
      % Create a new figure window
      delete(findall(0, '-depth',1, 'Tag','fpropertiesGUI'));
      hFig = figure('Number','off', 'Name','Application properties', 'Units','pixel', 'Pos',[300,200,500,500], 'Menu','none', 'Toolbar','none', 'Tag','fpropertiesGUI', 'Visible','off');
      hParent = hFig;
      setappdata(0,'isParamsGUIApproved',false)

      % Add the bottom action buttons
      btOK     = uicontrol('String','OK',     'Units','pixel', 'Pos',[ 50,5,60,30], 'Tag','btOK',     'Callback',@btOK_Callback);
      btCancel = uicontrol('String','Cancel', 'Units','pixel', 'Pos',[150,5,60,30], 'Tag','btCancel', 'Callback',@(h,e)close(hFig)); %#ok<NASGU>

      % Check the property values to determine whether the <OK> button should be enabled or not
      checkProps(propsList, btOK, true);
  
      % Set the figure icon & make visible
      jFrame = get(handle(hFig),'JavaFrame');
      icon = javax.swing.ImageIcon(fullfile(matlabroot, '/toolbox/matlab/icons/tool_legend.gif'));
      jFrame.setFigureIcon(icon);
      set(hFig, 'WindowStyle','modal', 'Visible','on');
      
      % Set the component's position
      %pos = [5,40,490,440];
      hFigPos = getpixelposition(hFig);
      pos = [5,40,hFigPos(3)-10,hFigPos(4)-50];
  else
      % Set the component's position
      drawnow;
      pos = getpixelposition(hParent);
      pos(1:2) = 5;
      pos = pos - [0,0,10,10];
      hFig = [];
  end

  %drawnow; pause(0.05);
  pane = javaObjectEDT(com.jidesoft.grid.PropertyPane(grid));
  customizePropertyPane(pane);
  [jPropsPane, hPropsPane_] = javacomponent(pane, pos, hParent);
  setappdata(hParent, 'jPropsPane',jPropsPane);
  setappdata(hParent, 'propsList',propsList);
  setappdata(hParent, 'propsHash',propsHash);
  setappdata(hParent, 'mirror',parameters);
  set(hPropsPane_,'tag','hpropertiesGUI');

  set(hPropsPane_, 'Units','norm');

  % Align the background colors
  bgcolor = pane.getBackground.getComponents([]);
  try set(hParent, 'Color', bgcolor(1:3)); catch, end  % this fails in uitabs - never mind (works ok in stand-alone figures)
  try pane.setBorderColor(pane.getBackground); catch, end  % error reported by Andrew Ness
  
  % If a new figure was created, make it modal and wait for user to close it
  if ~isempty(hFig)
      uiwait(hFig);
      if getappdata(0,'isParamsGUIApproved')
          parameters = test_data; %=getappdata(hFig, 'mirror');
      end
  end
  
  if nargout, hPropsPane = hPropsPane_; end  % prevent unintentional printouts to the command window
end  % propertiesGUI

%% Determine whether a specified object should be considered as having fields/properties
% Note: HG handles must be processed seperately for the main logic to work
function [hasProps,isHG] = hasProperties(object)
    % A bunch of tests, some of which may croak depending on the Matlab release, platform
    try isHG  = ishghandle(object); catch, isHG  = ishandle(object);  end
    try isst  = isstruct(object);   catch, isst  = false; end
    try isjav = isjava(object);     catch, isjav = false; end
    try isobj = isobject(object);   catch, isobj = false; end
    try isco  = iscom(object);      catch, isco  = false; end
    hasProps = ~isempty(object) && (isst || isjav || isobj || isco);
end

%% Customize the property-pane's appearance
function customizePropertyPane(pane)
  pane.setShowDescription(false);  % YMA: we don't currently have textual descriptions of the parameters, so no use showing an empty box that just takes up GUI space...
  pane.setShowToolBar(false);
  pane.setOrder(2);  % uncategorized, unsorted - see http://undocumentedmatlab.com/blog/advanced-jide-property-grids/#comment-42057
end

%% Prepare a list of some parameters for demo mode
function parameters = demoParameters
    parameters.floating_point_property = pi;
    parameters.signed_integer_property = int16(12);
    parameters.unsigned_integer_property = uint16(12);
    parameters.flag_property = true;
    parameters.file_property = mfilename('fullpath');
    parameters.folder_property = pwd;
    parameters.text_property = 'Sample text';
    parameters.fixed_choice_property = {'Yes','No','Maybe'};
    parameters.editable_choice_property = {'Yes','No','Maybe',''};  % editable if the last cell element is ''
    parameters.date_property = java.util.Date;  % today's date
    parameters.another_date_property = now-365;  % last year
    parameters.time_property = datestr(now,'HH:MM:SS');
    parameters.password_property = '*****';
    parameters.IP_address_property = '10.20.30.40';
    parameters.my_category.width = 4;
    parameters.my_category.height = 3;
    parameters.my_category.and_a_subcategory.is_OK = true;
    parameters.numeric_array_property = [11,12,13,14];
    parameters.cell_array_property  = {1,magic(3),'text',-4};
    parameters.color_property = [0.4,0.5,0.6];
    parameters.another_color_property = java.awt.Color.red;
    parameters.font_property = java.awt.Font('Arial', java.awt.Font.BOLD, 12);
    try parameters.class_object_property = matlab.desktop.editor.getActive; catch, end
end  % demoParameters

%% Prepare a list of properties
function propsList = preparePropsList(parameters, isEditable)
  propsList = java.util.ArrayList();

  % Convert a class object into a struct
  if isobject(parameters)
      parameters = struct(parameters);
  end

  % Check for an array of inputs (currently unsupported)
  %if numel(parameters) > 1,  error('YMA:propertiesGUI:ArrayParameters','Non-scalar inputs are currently unsupported');  end

  % Prepare a dynamic list of properties, based on the struct fields
  if isstruct(parameters) && ~isempty(parameters)
      %allParameters = parameters(:);  % convert ND array => 3D array
      allParameters = reshape(parameters, size(parameters,1),size(parameters,2),[]);
      numParameters = numel(allParameters);
      if numParameters > 1
          for zIdx = 1 : size(allParameters,3)
              for colIdx = 1 : size(allParameters,2)
                  for rowIdx = 1 : size(allParameters,1)
                      parameters = allParameters(rowIdx,colIdx,zIdx);
                      field_name = '';
                      field_label = sprintf('(%d,%d,%d)',rowIdx,colIdx,zIdx);
                      field_label = regexprep(field_label,',1\)',')');  % remove 3D if unnecesary
                      newProp = newProperty(parameters, field_name, field_label, isEditable, '', '', @propUpdatedCallback);
                      propsList.add(newProp);
                  end
              end
          end
      else
          % Dynamically (generically) inspect all the fields and assign corresponding props
          field_names = fieldnames(parameters);
          for field_idx = 1 : length(field_names)
              field_name = field_names{field_idx};
              value = parameters.(field_name);
              field_label = strrep(field_name,'_',' ');
              field_label(1) = upper(field_label(1));
              %if numParameters > 1,  field_label = [field_label '(' num2str(parametersIdx) ')'];  end
              field_description = '';  % TODO
              type = 'string';
              if isempty(value)
                  type = 'string';  % not really needed, but for consistency
              elseif isa(value,'java.awt.Color')
                  type = 'color';
              elseif isa(value,'java.awt.Font')
                  type = 'font';
              elseif isnumeric(value)
                  try %if length(value)==3
                      colorComponents = num2cell(value);
                      if numel(colorComponents) ~= 3
                          error(' ');  % bail out if definitely not a color
                      end
                      try
                          value = java.awt.Color(colorComponents{:});  % value between 0-1
                      catch
                          colorComponents = num2cell(value/255);
                          value = java.awt.Color(colorComponents{:});  % value between 0-255
                      end
                      type = 'color';
                  catch %else
                      if numel(value)==1
                          %value = value(1);
                          if value > now-3650 && value < now+3650
                              type = 'date';
                              value = java.util.Date(datestr(value));
                          elseif isa(value,'uint') || isa(value,'uint8') || isa(value,'uint16') || isa(value,'uint32') || isa(value,'uint64')
                              type = 'unsigned';
                          elseif isinteger(value)
                              type = 'signed';
                          else
                              type = 'float';
                          end
                      else
                          value = num2str(value);
                          if size(value,1) > size(value,2)
                              value = value';
                          end
                          if size(squeeze(value),2) > 1
                              % Convert multi-row string into a single-row string
                              value = [value'; repmat(' ',1,size(value,1))];
                              value = value(:)';
                          end
                          value = strtrim(regexprep(value,' +',' '));
                          if length(value) > 50
                              value(51:end) = '';
                              value = [value '...']; %#ok<AGROW>
                          end
                          value = ['[ ' value ' ]']; %#ok<AGROW>
                      end
                  end
              elseif islogical(value)
                  type = 'boolean';
              elseif ischar(value)
                  if exist(value,'dir')
                      type = 'folder';
                      value = java.io.File(value);
                  elseif exist(value,'file')
                      type = 'file';
                      value = java.io.File(value);
                  elseif value(1)=='*'
                      type = 'password';
                  elseif sum(value=='.')==3
                      type = 'IPAddress';
                  else
                      type = 'string';
                      if length(value) > 50
                          value(51:end) = '';
                          value = [value '...']; %#ok<AGROW>
                      end
                  end
              elseif iscellstr(value)
                  type = value;  % editable if the last cell element is ''
              elseif isa(value,'java.util.Date')
                  type = 'date';
              elseif isa(value,'java.io.File')
                  if value.isFile
                      type = 'file';
                  else  % value.isDirectory
                      type = 'folder';
                  end
              elseif iscell(value)
                  value = strtrim(regexprep(evalc('disp(value)'),' +',' '));
                  value = ['{ ' value ' }']; %#ok<AGROW>
              elseif isobject(value)
                  oldWarn = warning('off','MATLAB:structOnObject');
                  value = struct(value);
                  warning(oldWarn);
              elseif ~isstruct(value)
                  value = strtrim(regexprep(evalc('disp(value)'),' +',' '));
              end
              parameters.(field_name) = value;  % possibly updated above
              newProp = newProperty(parameters, field_name, field_label, isEditable, type, field_description, @propUpdatedCallback);
              propsList.add(newProp);
          end
      end
  else
      % You can also use direct assignments, instead of the generic code above. For example:
      % (Possible property types: signed, unsigned, float, file, folder, text or string, color, IPAddress, password, date, boolean, cell-array of strings)
      propsList.add(newProperty(parameters, 'flag_prop_name',   'Flag value:',     isEditable, 'boolean',            'Turn this on if you want to make extra plots', @propUpdatedCallback));
      propsList.add(newProperty(parameters, 'float_prop_name',  'Boolean prop',    isEditable, 'float',              'description 123...',   @propUpdatedCallback));
      propsList.add(newProperty(parameters, 'string_prop_name', 'My text msg:',    isEditable, 'string',             'Yaba daba doo',        @propUpdatedCallback));
      propsList.add(newProperty(parameters, 'int_prop_name',    'Now an integer',  isEditable, 'unsigned',           '123 456...',           @propUpdatedCallback));
      propsList.add(newProperty(parameters, 'choice_prop_name', 'And a drop-down', isEditable, {'Yes','No','Maybe'}, 'no description here!', @propUpdatedCallback));
  end
end  % preparePropsList

%% Prepare a data property
function prop = newProperty(dataStruct, propName, label, isEditable, dataType, description, propUpdatedCallback)

  % Auto-generate the label from the property name, if the label was not specified
  if isempty(label)
      label = strrep(propName,'_',' ');
      label(1) = upper(label(1));
  end

  % Create a new property with the chosen label
  prop = javaObjectEDT(com.jidesoft.grid.DefaultProperty);  % UNDOCUMENTED internal MATLAB component
  prop.setName(label);
  
  % Set the property to the current patient's data value
  try
      thisProp = dataStruct.(propName);
  catch
      thisProp = dataStruct;
  end
  origProp = thisProp;
  if isstruct(thisProp)  %hasProperties(thisProp)
      % Accept any object having data fields/properties
      try
          thisProp = get(thisProp);
      catch
          oldWarn = warning('off','MATLAB:structOnObject');
          thisProp = struct(thisProp);
          warning(oldWarn);
      end

      % Parse the children props and add them to this property
      %summary = regexprep(evalc('disp(thisProp)'),' +',' ');
      %prop.setValue(summary);  % TODO: display summary dynamically
      if numel(thisProp) < 2
          prop.setValue('');
      else
          sz = size(thisProp);
          szStr = regexprep(num2str(sz),' +','x');
          prop.setValue(['[' szStr ' struct array]']);
      end
      prop.setEditable(false);
      children = toArray(preparePropsList(thisProp, isEditable));
      for childIdx = 1 : length(children)
          prop.addChild(children(childIdx));
      end
  else
      prop.setValue(thisProp);
      prop.setEditable(isEditable);
  end

  % Set property editor, renderer and alignment
  if iscell(dataType)
      % treat this as drop-down values
      cbIsEditable = true;
      if isempty(dataType{end})  % ends with '' - editable
          dataType(end) = [];  % remove from the drop-down list
      else  % standard drop-down, non-editable
          cbIsEditable = false;
      end
      editor = com.jidesoft.grid.ListComboBoxCellEditor(dataType);
      try editor.getComboBox.setEditable(cbIsEditable); catch, end % #ok<NOCOM>
      %set(editor,'EditingStoppedCallback',{@propUpdatedCallback,tagName,propName});
      alignProp(prop, editor);
      try prop.setValue(origProp{1}); catch, end
  else
      switch lower(dataType)
          case 'signed',    %alignProp(prop, com.jidesoft.grid.IntegerCellEditor,    'int32');
                            model = javax.swing.SpinnerNumberModel(prop.getValue, -intmax, intmax, 1);
                            editor = com.jidesoft.grid.SpinnerCellEditor(model);
                            alignProp(prop, editor, 'int32');
          case 'unsigned',  %alignProp(prop, com.jidesoft.grid.IntegerCellEditor,    'uint32');
                            val = max(0, min(prop.getValue, intmax));
                            model = javax.swing.SpinnerNumberModel(val, 0, intmax, 1);
                            editor = com.jidesoft.grid.SpinnerCellEditor(model);
                            alignProp(prop, editor, 'uint32');
          case 'float',     alignProp(prop, com.jidesoft.grid.CalculatorCellEditor, 'double');  % DoubleCellEditor
          case 'boolean',   alignProp(prop, com.jidesoft.grid.BooleanCheckBoxCellEditor, 'logical');
          case 'folder',    alignProp(prop, com.jidesoft.grid.FolderCellEditor);
          case 'file',      alignProp(prop, com.jidesoft.grid.FileCellEditor);
          case 'ipaddress', alignProp(prop, com.jidesoft.grid.IPAddressCellEditor);
          case 'password',  alignProp(prop, com.jidesoft.grid.PasswordCellEditor);
          case 'color',     alignProp(prop, com.jidesoft.grid.ColorCellEditor);
          case 'font',      alignProp(prop, com.jidesoft.grid.FontCellEditor);
          case 'text',      alignProp(prop);
          case 'time',      alignProp(prop);  % maybe use com.jidesoft.grid.FormattedTextFieldCellEditor ?

          case 'date',      dateModel = com.jidesoft.combobox.DefaultDateModel;
                            dateFormat = java.text.SimpleDateFormat('dd/MM/yyyy');
                            dateModel.setDateFormat(dateFormat);
                            editor = com.jidesoft.grid.DateCellEditor(dateModel, 1);
                            alignProp(prop, editor, 'java.util.Date');
                            try
                                prop.setValue(dateFormat.parse(prop.getValue));  % convert string => Date
                            catch
                                % ignore
                            end

          otherwise,        alignProp(prop);  % treat as a simple text field
      end
  end  % for all possible data types

  prop.setDescription(description);
  if ~isempty(description)
      renderer = com.jidesoft.grid.CellRendererManager.getRenderer(prop.getType, prop.getEditorContext);
      renderer.setToolTipText(description);
  end

  % Set the property's editability state
  if prop.isEditable
      % Set the property's label to be black
      prop.setDisplayName(['<html><font size="4" color="black">' label]);

      % Add callbacks for property-change events
      hprop = handle(prop, 'CallbackProperties');
      set(hprop,'PropertyChangeCallback',{propUpdatedCallback,propName});
  else
      % Set the property's label to be gray
      prop.setDisplayName(['<html><font size="4" color="gray">' label]);
  end
  
  setPropName(prop,propName);
end  % newProperty

%% Set property name in the Java property reference
function setPropName(hProp,propName)
    try
        set(hProp,'UserData',propName)
    catch
        %setappdata(hProp,'UserData',propName)
        hp = schema.prop(handle(hProp),'UserData','mxArray'); %#ok<NASGU>
        set(handle(hProp),'UserData',propName)
    end
end  % setPropName

%% Get property name from the Java property reference
function propName = getPropName(hProp)
    try
        propName = get(hProp,'UserData');
    catch
        %propName = char(getappdata(hProp,'UserData'));
        propName = get(handle(hProp),'UserData');
    end
end  % getPropName

%% Align a text property to right/left
function alignProp(prop, editor, propTypeStr, direction)
  if nargin < 2 || isempty(editor),      editor = com.jidesoft.grid.StringCellEditor;  end  %(javaclass('char',1));
  if nargin < 3 || isempty(propTypeStr), propTypeStr = 'cellstr';  end  % => javaclass('char',1)
  if nargin < 4 || isempty(direction),   direction = javax.swing.SwingConstants.RIGHT;  end

  % Set this property's data type
  propType = javaclass(propTypeStr);
  prop.setType(propType);

  % Prepare a specific context object for this property
  if strcmpi(propTypeStr,'logical')
      %TODO - FIXME
      context = editor.CONTEXT;
      prop.setEditorContext(context);
      %renderer = CheckBoxRenderer;
      %renderer.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
      %com.jidesoft.grid.CellRendererManager.registerRenderer(propType, renderer, context);
  else
      context = com.jidesoft.grid.EditorContext(prop.getName);
      prop.setEditorContext(context);

      % Register a unique cell renderer so that each property can be modified seperately
      %renderer = com.jidesoft.grid.CellRendererManager.getRenderer(propType, prop.getEditorContext);
      renderer = com.jidesoft.grid.ContextSensitiveCellRenderer;
      com.jidesoft.grid.CellRendererManager.registerRenderer(propType, renderer, context);
      renderer.setBackground(java.awt.Color.white);
      renderer.setHorizontalAlignment(direction);
      %renderer.setHorizontalTextPosition(direction);
  end

  % Update the property's cell editor
  try editor.setHorizontalAlignment(direction); catch, end
  try editor.getTextField.setHorizontalAlignment(direction); catch, end
  try editor.getComboBox.setHorizontalAlignment(direction); catch, end
  
  % Set limits on unsigned int values
  try
      if strcmpi(propTypeStr,'uint32')
          %pause(0.01);
          editor.setMinInclusive(java.lang.Integer(0));
          editor.setMinExclusive(java.lang.Integer(-1));
          editor.setMaxExclusive(java.lang.Integer(intmax));
          editor.setMaxInclusive(java.lang.Integer(intmax));
      end
  catch
      % ignore
  end
  com.jidesoft.grid.CellEditorManager.registerEditor(propType, editor, context);
end  % alignProp

%% Property updated callback function
function propUpdatedCallback(prop, eventData, propName)
  try if strcmpi(char(eventData.getPropertyName),'parent'),  return;  end;  catch, end
  hFig = findall(0, '-depth',1, 'Tag','fpropertiesGUI');
  if isempty(hFig)
      hPropsPane = findall(0,'Tag','hpropertiesGUI');
      if isempty(hPropsPane),  return;  end
      hFig = get(hPropsPane,'Parent');
  end
  if isempty(hFig),  return;  end
  propsList = getappdata(hFig, 'propsList');
  propsPane = getappdata(hFig, 'jPropsPane');
  data = getappdata(hFig, 'mirror');

  % Get the updated property value
  propValue = get(prop,'Value');
  if isjava(propValue)
      if isa(propValue,'java.util.Date')
          sdf = java.text.SimpleDateFormat('MM-dd-yyyy');
          propValue = datenum(sdf.format(propValue).char);  %#ok<NASGU>
      elseif isa(propValue,'java.awt.Color')
          propValue = propValue.getColorComponents([])';  %#ok<NASGU>
      else
          propValue = char(propValue);  %#ok<NASGU>
      end
  end

  % Get the actual recursive propName
  try
      oldWarn = warning('off','MATLAB:hg:JavaSetHGProperty');
      try prop = java(prop); catch, end
      while isa(prop,'com.jidesoft.grid.Property')
          prop = get(prop,'Parent');
          newName = getPropName(prop);
          if isempty(newName), break; end
          propName = [newName '.' propName]; %#ok<AGROW>
      end
  catch
      % Reached the top of the property's heirarchy - bail out
      warning(oldWarn);
  end
  
  % Update the mirror with the updated field value
  %data.(propName) = propValue;  % croaks on multiple sub-fields
  eval(['data.' propName ' = propValue;']);

  % Update the local mirror
  setappdata(hFig, 'mirror',data);
  
  % Update the display
  checkProps(propsList, hFig);
  try propsPane.repaint; catch; end
end  % propUpdatedCallback

%% <OK> button callback function
function btOK_Callback(btOK, eventData) %#ok<INUSD>
  global test_data

  % Store the current data-info struct mirror in the global struct
  hFig = ancestor(btOK, 'figure');
  test_data = getappdata(hFig, 'mirror');
  setappdata(0,'isParamsGUIApproved',true);

  % Close the window
  try
      close(hFig);
  catch
      delete(hFig);  % force-close
  end
end  % btOK_Callback

%% Check whether all mandatory fields have been filled, update background color accordingly
function checkProps(propsList, hContainer, isInit)
    if nargin < 3,  isInit = false;  end
    okEnabled = 'on';
    try propsArray = propsList.toArray(); catch, return; end
    for propsIdx = 1 : length(propsArray)
        isOk = checkProp(propsArray(propsIdx));
        if ~isOk || isInit,  okEnabled = 'off';  end
    end
    
    % Update the <OK> button's editability state accordingly
    btOK = findall(hContainer, 'Tag','btOK');
    set(btOK, 'Enable',okEnabled);
    drawnow; pause(0.01);
end  % checkProps

function isOk = checkProp(prop)
  isOk = true;
  oldWarn = warning('off','MATLAB:hg:JavaSetHGProperty');
  warning off MATLAB:hg:PossibleDeprecatedJavaSetHGProperty
  propName = getPropName(prop);
  renderer = com.jidesoft.grid.CellRendererManager.getRenderer(get(prop,'Type'), get(prop,'EditorContext'));
  warning(oldWarn);
  mandatoryFields = {};  % TODO - add the mandatory field-names here
  if any(strcmpi(propName, mandatoryFields)) && isempty(get(prop,'Value'))
      propColor = java.awt.Color.yellow;
      isOk = false;
  elseif ~prop.isEditable
      %propColor = java.awt.Color.gray;
      propColor = renderer.getBackground();
  else
      propColor = java.awt.Color.white;
  end
  renderer.setBackground(propColor);
end  % checkProp

%% Return java.lang.Class instance corresponding to the Matlab type
function jclass = javaclass(mtype, ndims)
    % Input arguments:
    % mtype:
    %    the MatLab name of the type for which to return the java.lang.Class
    %    instance
    % ndims:
    %    the number of dimensions of the MatLab data type
    %
    % See also: class
    
    % Copyright 2009-2010 Levente Hunyadi
    % Downloaded from: http://www.UndocumentedMatlab.com/files/javaclass.m
    
    validateattributes(mtype, {'char'}, {'nonempty','row'});
    if nargin < 2
        ndims = 0;
    else
        validateattributes(ndims, {'numeric'}, {'nonnegative','integer','scalar'});
    end
    
    if ndims == 1 && strcmp(mtype, 'char');  % a character vector converts into a string
        jclassname = 'java.lang.String';
    elseif ndims > 0
        jclassname = javaarrayclass(mtype, ndims);
    else
        % The static property .class applied to a Java type returns a string in
        % MatLab rather than an instance of java.lang.Class. For this reason,
        % use a string and java.lang.Class.forName to instantiate a
        % java.lang.Class object; the syntax java.lang.Boolean.class will not do so
        switch mtype
            case 'logical'  % logical vaule (true or false)
                jclassname = 'java.lang.Boolean';
            case 'char'  % a singe character
                jclassname = 'java.lang.Character';
            case {'int8','uint8'}  % 8-bit signed and unsigned integer
                jclassname = 'java.lang.Byte';
            case {'int16','uint16'}  % 16-bit signed and unsigned integer
                jclassname = 'java.lang.Short';
            case {'int32','uint32'}  % 32-bit signed and unsigned integer
                jclassname = 'java.lang.Integer';
            case {'int64','uint64'}  % 64-bit signed and unsigned integer
                jclassname = 'java.lang.Long';
            case 'single'  % single-precision floating-point number
                jclassname = 'java.lang.Float';
            case 'double'  % double-precision floating-point number
                jclassname = 'java.lang.Double';
            case 'cellstr'  % a single cell or a character array
                jclassname = 'java.lang.String';
            otherwise
                jclassname = mtype;
                %error('java:javaclass:InvalidArgumentValue', ...
                %    'MatLab type "%s" is not recognized or supported in Java.', mtype);
        end
    end
    % Note: When querying a java.lang.Class object by name with the method
    % jclass = java.lang.Class.forName(jclassname);
    % MatLab generates an error. For the Class.forName method to work, MatLab
    % requires class loader to be specified explicitly.
    jclass = java.lang.Class.forName(jclassname, true, java.lang.Thread.currentThread().getContextClassLoader());
end  % javaclass
    
%% Return the type qualifier for a multidimensional Java array
function jclassname = javaarrayclass(mtype, ndims)
    switch mtype
        case 'logical'  % logical array of true and false values
            jclassid = 'Z';
        case 'char'  % character array
            jclassid = 'C';
        case {'int8','uint8'}  % 8-bit signed and unsigned integer array
            jclassid = 'B';
        case {'int16','uint16'}  % 16-bit signed and unsigned integer array
            jclassid = 'S';
        case {'int32','uint32'}  % 32-bit signed and unsigned integer array
            jclassid = 'I';
        case {'int64','uint64'}  % 64-bit signed and unsigned integer array
            jclassid = 'J';
        case 'single'  % single-precision floating-point number array
            jclassid = 'F';
        case 'double'  % double-precision floating-point number array
            jclassid = 'D';
        case 'cellstr'  % cell array of strings
            jclassid = 'Ljava.lang.String;';
        otherwise
            jclassid = ['L' mtype ';'];
            %error('java:javaclass:InvalidArgumentValue', ...
            %    'MatLab type "%s" is not recognized or supported in Java.', mtype);
    end
    jclassname = [repmat('[',1,ndims), jclassid];
end  % javaarrayclass
