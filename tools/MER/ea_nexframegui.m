function varargout = ea_nexframegui(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ea_nexframegui_OpeningFcn, ...
                   'gui_OutputFcn',  @ea_nexframegui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ea_nexframegui is made visible.
function ea_nexframegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ea_nexframegui (see VARARGIN)

% Choose default command line output for ea_nexframegui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

merframe = varargin{1};
setappdata(hObject, 'merframe', merframe);

gui_build_table(handles);

set(handles.figure1,'WindowStyle','modal')
uiwait(handles.figure1);

% UIWAIT makes ea_nexframegui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function gui_build_table(handles)

side_strings = {'right', 'left'};
marker_strings = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'Entry'};
coord_strings = {'X', 'Y', 'Z'};

merframe = getappdata(handles.figure1, 'merframe');

% Destroy existing table
n_rows = count_table_rows(handles);
if n_rows > 0
    ui_tags = {handles.figure1.Children.Tag};
    for row_ix = 1:n_rows
        tag = ['popup_side_' num2str(row_ix)];
        delete(handles.figure1.Children(strcmpi(ui_tags, tag)));
        tag = ['popup_marker_' num2str(row_ix)];
        delete(handles.figure1.Children(strcmpi(ui_tags, tag)));
        for dim_id = 1:3
            tag = ['edit_pos_' coord_strings{dim_id} '_' num2str(row_ix)];
            delete(handles.figure1.Children(strcmpi(ui_tags, tag)));
        end
    end
    clear ui_tags row_ix tag dim_id
end

% Build new table
side_x = handles.colheader_hemisphere.Position(1);
top_y = handles.colheader_hemisphere.Position(2);
marker_x = handles.colheader_landmark.Position(1);
dim_x = [   handles.colheader_X.Position(1),...
            handles.colheader_Y.Position(1),...
            handles.colheader_Z.Position(1)];
col_width = 70;
row_height = 20;


row_ix = 0;
for side_ix = 1:length(merframe)
    sid = length(merframe) + 1 - side_ix;  % reverse order
    markers = merframe{sid}.markers;
    for marker_ix = 1:length(merframe{sid}.markers)
        row_ix = row_ix + 1;
        side_vals(row_ix) = sid;
        marker_vals(row_ix) = find(strcmpi(marker_strings, merframe{sid}.markers{marker_ix}));
        dim_vals(row_ix, :) = merframe{sid}.coords(marker_ix, :);
    end
end
    
for row_ix = 1:length(side_vals)
    row_y = top_y - 30*row_ix;
    % Add popup for side
    uicontrol(handles.figure1,...
        'Style', 'popupmenu',...
        'Tag', ['popup_side_' num2str(row_ix)],...
        'String', side_strings,...
        'Value', side_vals(row_ix),...
        'Position', [side_x row_y col_width row_height],...
        'Units', 'pixels',...
        'BackgroundColor', [1 1 1],...
        'Enable', 'off');  % TODO: Enable when support for other configs is added.
    % Add popup for marker A-H, Entry
    uicontrol(handles.figure1,...
        'Style', 'popupmenu',...
        'Tag', ['popup_marker_' num2str(row_ix)],...
        'String', marker_strings,...
        'Value', marker_vals(row_ix),...
        'Position', [marker_x row_y col_width row_height],...
        'Units', 'pixels',...
        'BackgroundColor', [1 1 1],...
        'Enable', 'off');  % TODO: Enable when support for other configs is added.
    for dim_id = 1:3
        % Add edit box for X, Y, Z
        uicontrol(handles.figure1,...
            'Style', 'edit',...
            'Tag', ['edit_pos_' coord_strings{dim_id} '_' num2str(row_ix)],...
            'String', num2str(dim_vals(row_ix, dim_id)),...
            'Position', [dim_x(dim_id) row_y col_width row_height],...
            'Units', 'pixels',...
            'BackgroundColor', [1 1 1]);
    end
end


% --- Outputs from this function are returned to the command line.
function varargout = ea_nexframegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = getappdata(handles.figure1, 'merframe');
delete(handles.figure1);


% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pb_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% - take the norm of each vector to get frame base coordinate system
merframe = getappdata(handles.figure1, 'merframe');

ui_tags = {handles.figure1.Children.Tag};
n_rows = count_table_rows(handles);
dim_str = {'X', 'Y', 'Z'};
sides = nan(n_rows, 1);
markers = cell(n_rows, 1);
coords = nan(n_rows, 3);
for row_ix = 1:n_rows
    side_tag = ['popup_side_' num2str(row_ix)];
    hside = handles.figure1.Children(strcmpi(ui_tags, side_tag));
    sides(row_ix) = hside.Value;
    marker_tag = ['popup_marker_' num2str(row_ix)];
    hmarker = handles.figure1.Children(strcmpi(ui_tags, marker_tag));
    markers{row_ix} = hmarker.String{hmarker.Value};
    for dim_ix = 1:3
        dim_tag = ['edit_pos_' dim_str{dim_ix} '_' num2str(row_ix)];
        hdim = handles.figure1.Children(strcmpi(ui_tags, dim_tag));
        coords(row_ix, dim_ix) = str2double(hdim.String);
    end
end

% Store in merstruct
for side_ix = 1:2
    b_rows = sides == side_ix;
    merframe{side_ix}.markers = markers(b_rows)';
    merframe{side_ix}.coords = coords(b_rows, :);
end

% TODO: Support other configurations than requiring A, E, and Entry.
AEP_str = {'A', 'E', 'Entry'};
for side_ix = 1:2
    AEP = nan(3, 3);
    for aep_ix = 1:3
        b_row = strcmpi(markers, AEP_str{aep_ix}) & sides == side_ix;
        AEP(aep_ix, :) = coords(b_row, :);
    end
    m = (AEP(1, :) + AEP(2, :)) / 2;  %midpoint between A and E becomes origin.
    y = AEP(1, :) - m;  % Vector from frame origin to frame +y in scrf space
    z = m - AEP(3, :);  % '' for +z
    x = cross(y, z);    % '' for +x
    x = x / norm(x); y = y / norm(y); z = z / norm(z);  % Make unit vectors.
    
    % Get the transformation between coordinate systems
    def_spc = [1 0 0; 0 1 0; 0 0 1];
    xform = [x; y; z] \ def_spc;  % == inv([x; y; z])
    % xform gives full nexframe transform in scrf space.
    % However, lead localization constrains pitch and roll, so we only need
    % yaw.
    % Save yaw in merframe
    merframe{side_ix}.yaw_rad = atan2(xform(2, 1), xform(1, 1));
end

setappdata(handles.figure1, 'merframe', merframe);
uiresume(handles.figure1);

% TODO: I'm not positive yaw-only is correct. We might have to take
% our original xform, transform to mni space, then do a best fit alignment
% of its z-vec with the DBS lead localization.

% % Visualization:
% i_ = [1 0 0]; j_ = [0 1 0]; k_ = [0 0 1]; n_ = [0 0 0]; s_ = 3;
% quiver3(n_, n_, n_, s_*i_, s_*j_, s_*k_, 'k');
% hold on
% quiver3(n_, n_, n_, s_*x, s_*y, s_*z, 'LineWidth', 3);
% xlim([-s_, s_]); ylim([-s_, s_]); zlim([-s_, s_]);
% xlabel('x'); ylabel('y'); zlabel('z');
% scatter3(mer_offsets(:,1), mer_offsets(:,2), mer_offsets(:,3), 'k');
% scatter3(full_rot_offsets(:,1), full_rot_offsets(:,2), full_rot_offsets(:,3), 'r');
% scatter3(yaw_only_offsets(:,1), yaw_only_offsets(:,2), yaw_only_offsets(:,3), 'k', 'filled');



% --- Executes on button press in pb_cancel.
function pb_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end


function n_rows = count_table_rows(handles)
ui_tags = {handles.figure1.Children.Tag};
n_rows = 0;
for tag_ix = 1:length(ui_tags)
    if (length(ui_tags{tag_ix}) > 11) &&...
        strcmpi(ui_tags{tag_ix}(1:11), 'popup_side_')
        n_rows = max(n_rows, str2double(handles.figure1.Children(tag_ix).Tag(12)));
    end
end
