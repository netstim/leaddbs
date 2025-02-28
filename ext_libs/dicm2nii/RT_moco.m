function RT_moco()
% Display and save motion information at real time. It also shows images and the
% progress of scanning, and allows to check motion information for previous
% series/patients.
%
% To make this work, you will need:
%  1. Set up shared folder at the computer running RT_moco. 
%     The folder default to ../incoming_DICOM, but better set to your own by 
%     setpref('dicm2nii_gui_para', 'incomingDcm', '/mypath/myIncomingDicom');
%     The result log will be saved into incoming_DICOM/RTMM_log/ folder.
%  2. Set up real time image transfer at Siemens console.

% 200207 xiangrui.li at gmail.com first working version inspired by FIRMM
% 240115 works for XA30 and drop support for mosaic

% Create/re-use GUI and start timer
fh = findall(0, 'Type', 'figure', 'Tag', 'RT_moco');
if ~isempty(fh), figure(fh); return; end
res = get(0, 'ScreenSize');
fh = figure('mc'*[256 1]'); clf(fh);
set(fh, 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', ... 
	'DockControls', 'off', 'CloseRequestFcn', @closeFig, 'Color', [1 1 1]*0.94, ...
    'Name', 'Real Time Image Monitor ', 'Tag', 'RT_moco', ...
    'WindowState', 'maximized', 'Visible', 'off');
hs.fig = fh;
hs.rootDir = getpref('dicm2nii_gui_para', 'incomingDcm', '../incoming_DICOM/');
hs.backupDir = getpref('dicm2nii_gui_para', 'backupDir', '');
fullName = dicm2nii('', 'fullName', 'func_handle');
hs.rootDir = fullName(hs.rootDir);
if isfolder(hs.backupDir), hs.backupDir = fullName(hs.backupDir); end
hs.logDir = [hs.rootDir 'RTMM_log/'];
if ~isfolder(hs.logDir), mkdir(hs.logDir); end % folder to save subj.mat

h = uimenu(fh, 'Label', '&Patient');
hs.menu(1) = uimenu(h, 'Label', 'Load Patient', 'Callback', @loadSubj);
hs.menu(2) = uimenu(h, 'Label', 'Redo Patient', 'Callback', @redoSubj);
hs.menu(3) = uimenu(h, 'Label', 'Create QC Report', 'Callback', @QC_report);
hs.menu(4) = uimenu(h, 'Label', 'Close Patient', 'Callback', @closeSubj);

h = uimenu(fh, 'Label', '&Series');
uimenu(h, 'Label', 'View Selected Series in 3D', 'Callback', @view_3D);
uimenu(h, 'Label', 'Overlay Selected Series onto Anatomy', 'Callback', @overlay);
hs.derived = uimenu(h, 'Label', 'Skip DERIVED Series', 'Checked', 'on', ...
    'Callback', @toggleChecked, 'Separator', 'on');
hs.SBRef = uimenu(h, 'Label', 'Skip *_SBRef Series', 'Callback', @toggleChecked, 'Checked', 'on');

h = uimenu(fh, 'Label', '&View');
uimenu(h, 'Label', 'Reset Brightness', 'Callback', @setCLim);
uimenu(h, 'Label', 'Increase Brightness', 'Callback', @setCLim);
uimenu(h, 'Label', 'Decrease Brightness', 'Callback', @setCLim);
hDV = uimenu(h, 'Label', '&DVARS Threshold', 'Separator', 'on');
for i = [0.01 0.02 0.04 0.05 0.06 0.08 0.1 0.2 0.4]
    uimenu(hDV, 'Label', num2str(i), 'Callback', @DV_yLim);
end
uimenu(h, 'Label', 'Show FD plot', 'Callback', @toggleFD, 'Separator', 'on');
hFD = uimenu(h, 'Label', '&FD Axis Range');
for i = [0.18 0.3:0.3:1.5 2.4 3 6]
    uimenu(hFD, 'Label', num2str(i), 'Callback', @FD_range)
end

panel = @(pos)uipanel(fh, 'Position', pos, 'BorderType', 'none');
if res(3) < res(4) % Portrait
    pa1 = panel([0 0.62 1 0.38]); % img and label
    pa2 = panel([0 0 1 0.62]); % table and plot
    axPos = [0.05 0 0.65 1]; % img axis
    lbPos = [0.65 0 0.34 1]; % label axis
    subjPos = [1 0.95]; seriesPos = [1 0.1]; msPos = [1 0.7]; ha = 'right';
else
    pa1 = panel([0 0 0.38 1]);
    pa2 = panel([0.38 0 0.62 1]);
    axPos = [0.05 0.31 0.9 0.67];
    lbPos = [0.05 0.01 0.9 0.3]; 
    subjPos = [0 0.98]; seriesPos = [1 0.1]; msPos = [0 0.6]; ha = 'left';
end

dy = 0.06 * (0:3);
hs.ax = axes(pa2, 'Position', [0.07 0.5 0.86 0.38], ...
    'NextPlot', 'add', 'XLim', [0.5 300.5], 'UserData', dy, ...
    'TickDir', 'out', 'TickLength', 0.002*[1 1], 'ColorOrder', [0 0 1; 1 0 1]);
xlabel(hs.ax, 'Volume Number');
hs.slider = uicontrol(pa2, 'Units', 'normalized', 'Position', [0.05 0.96 0.9 0.03], ...
    'Style', 'slider', 'Value', 1, 'Min', 1, 'Max', 300, 'Callback', @sliderCB, ...
    'BackgroundColor', 0.5*[1 1 1], 'SliderStep', [1 1]./300);

yyaxis left; ylabel(hs.ax, 'DVARS');
set(hs.ax, 'YTick', dy, 'YLim', dy([1 4]));
c3 = [0 0.8 0;  0.8 0.8 0;  0.5 0 0];
for i = 3:-1:1
    rectangle(hs.ax, 'Position', [0.5 dy(i) 2000 dy(i+1)-dy(i)], ...
        'FaceColor', c3(i,:), 'EdgeColor', c3(i,:), 'LineWidth', 0.01);
end
hs.dv = plot(hs.ax, nan, '.:');

yyaxis right; ylabel(hs.ax, 'Framewise Displacement (mm)');
set(hs.ax, 'YTick', 0:0.4:1.2, 'YLim', [0 1.2]);

txt = @(a)text(hs.ax, 'Units', 'normalized', 'Position', a, 'FontSize', 12, ...
    'BackgroundColor', [1 1 1]*0.94, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
hs.pct(1) = txt([0.995 0.32]); hs.pct(2) = txt([0.995 0.65]);

hs.fd = plot(hs.ax, 0, '.:', 'Visible', 'off');
hs.ax.YAxis(2).Visible = 'off';

ax = axes(pa2, 'Position', [0.07 0.88 0.86 0.04], 'Color', fh.Color);
hs.resp = plot(ax, nan, ones(1,3), 'o', 'Color', 'none', 'MarkerSize', 8);
set(hs.resp(1), 'MarkerFaceColor', [0 0 0]);
set(hs.resp(2), 'MarkerFaceColor', [1 0 0]);
set(hs.resp(3), 'MarkerFaceColor', [0 0.8 0]);
title(ax, ' ', 'FontSize', 14, 'interpreter', 'tex');
set(ax, 'XLim', [0.5 300.5], 'YLim', [0.5 1.5], 'Visible', 'off');

vars = {'Description' 'Series' 'Instances' '<font color="#00cc00">Green</font>' ...
    '<font color="#cccc00">Yellow</font>' 'MeanFD'};
w2 = [90 100 90 90 100]; w2 = num2cell([res(3)*pa2.Position(3)*0.94-sum(w2)-24 w2]);
hs.table = uitable(pa2, 'Units', 'normalized', 'Position', [0.02 0.01 0.96 0.42], ...
    'FontSize', 14, 'RowName', [], 'CellSelectionCallback', @tableCB, ...
    'ColumnName', strcat('<html><h2>', vars, '</h2></html>'), 'ColumnWidth', w2);

ax = axes(pa1, 'Position', axPos, 'YDir', 'reverse', 'Visible', 'off', 'CLim', [0 1]);
hs.img = image(ax, 'CData', ones(2)*0.94, 'CDataMapping', 'scaled');
axis equal; colormap(ax, 'gray');
hs.instnc = text(ax, 'Units', 'normalized', 'Position', [0.99 0.01], ...
    'Color', 'y', 'FontSize', 14, 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'bottom');

ax = axes(pa1, 'Position', lbPos, 'Visible', 'off');
hs.subj = text(ax, 'Position', subjPos, 'FontSize', 24, 'FontWeight', 'bold', ...
    'HorizontalAlignment', ha, 'VerticalAlignment', 'top', 'Interpreter', 'none');
hs.series = text(ax, 'Position', seriesPos, 'FontSize', 18, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
    'BackgroundColor', fh.Color, 'Interpreter', 'none');
hs.MMSS = text(ax, 'Position', msPos, 'FontSize', 18, 'FontWeight', 'bold', ...
    'HorizontalAlignment', ha, 'VerticalAlignment', 'top', 'Interpreter', 'none', ...
    'BackgroundColor', fh.Color, 'String', {'' ''}, 'Color', 'b');

set(fh, 'HandleVisibility', 'off', 'Visible', 'on'); % fh.Resize = 'off';

% set up serial port
if ispc, port = 'COM1'; else, port = '/dev/ttyUSB0'; end % change this to yours
try %#ok<*TRYNC>
    hs.serial = serialport(port, 115200, 'Timeout', 0.3);
    hs.serial.UserData = struct('fig', fh, 'send', false);
    configureCallback(hs.serial, "byte", 1, @serialRead);  
end

hs.timer = timer('StartDelay', 5, 'ObjectVisibility', 'off', 'UserData', fh, ...
    'StopFcn', @saveResult, 'TimerFcn', @timerFcn, 'ErrorFcn', @errorLog);

hs.countDown = timer('ExecutionMode', 'fixedRate', 'ObjectVisibility', 'off', ...
    'TimerFcn', {@countDown hs}, 'UserData', datetime);

guidata(fh, hs); closeSubj(fh);
start(hs.timer);

%% Log error for debugging
function errorLog(obj, evt) %#ok
hs = guidata(obj.UserData);
vnam = hs.series.String;
if iscell(vnam) && numel(vnam)>1, vnam = vnam{1}; end
vnam = hs.subj.String+"_"+hs.series.UserData+hs.instnc.String+"_"+vnam;
vnam = genvarname("err_"+vnam);
eval(vnam +" = evt;");
fnam = hs.logDir+"errorLog.mat";
if isfile(fnam), save(fnam, vnam, '-append'); else, save(fnam, vnam); end

%% TimerFunc: do a series if avail, then call stopFunc to save result
function timerFcn(obj, ~)
try doSeries(obj); catch me, errorLog(obj, me); end
function doSeries(obj)
hs = guidata(obj.UserData);
hs.fig.Name(end) = 78 - hs.fig.Name(end); % dot/space switch: timer indicator
if hs.timer.StartDelay > 1, return; end % no new series
set(hs.menu, 'Enable', 'off'); set([hs.table hs.slider], 'Enable', 'inactive');

bNam = fullfile(hs.subj.UserData, hs.series.UserData);
nam = dir([bNam '000001_*.dcm']);
if isempty(nam), return; end
s = dicm_hdr_wait(nam, []);
if isempty(s), return; end % non-image dicom, skip series
if size(hs.table.Data,1)<1 && datetime-datetime(nam.date)<minutes(5)
    save([hs.rootDir '/hdr_' hs.subj.String], 's'); % as new subj flag
end

hs.fig.WindowState = 'maximized';
L0 = java.awt.MouseInfo.getPointerInfo().getLocation();
java.awt.Robot().mouseMove(L0.getX+9, L0.getY+9); % wake up screen
pause(0.1); java.awt.Robot().mouseMove(L0.getX, L0.getY);

if hs.derived.Checked=="on" && contains(s.ImageType, 'DERIVED'), return; end
if hs.SBRef.Checked=="on" && endsWith(s.SeriesDescription, '_SBRef'), return; end
if startsWith(s.SequenceName, 'ABCD3d1'), return; end
if strcmp(s.SeriesDescription, 'MoCoSeries'), return; end

nTE = asc_header(s, 'lContrasts', 1);
% fieldmap phase diff series: nTE=2, EchoNumber=2
if (contains(s.ImageType, '\P\') && contains(s.SequenceName, 'fm2d')) || ...
   (contains(s.ImageType, '\MEAN') && s.NumberOfAverages>1) % vNAV RMS
    nTE = 1;
end

multiFrameFields = dicm2nii('', 'multiFrameFields', 'func_handle');
s = multiFrameFields(s);
nCh = 1;
if s.ICE_Dims(1) ~= "X" % coil not combined: no \CC: in ImageHistory
    k = 'sCoilSelectMeas.aRxCoilSelectData[0].asList.__attribute__.size';
    nCh = asc_header(s, k, 1);
end

isDTI = contains(s.ImageType, '\DIFFUSION');
if isDTI
    nTR = [];
    for i = 1:999
        a = asc_header(s, sprintf('sDiffusion.alAverages[%i]', i-1));
        if isempty(a), break; else, nTR(i) = a; end %#ok
    end
    nTR = nTR(1) + sum(nTR(2:end))*asc_header(s, 'sDiffusion.lDiffDirections');
else % ToDo: Reduce to real nTR for Sparse seq
    nTR = asc_header(s, 'lRepetitions', 0) + 1;
end
hs.series.String = seriesInfo(s);
init_series(hs, s, nTR);
img = dicm_img(s);
set_img(hs.img, img);
drawnow; 

% T1, T2, fieldmap etc: show info/img only
if ~isDTI && ~contains(s.SequenceName, {'epfid2d' 'ASL'}) || nCh>1
    hs.instnc.String = '1';
    hs.slider.Value = 1;
    return;
end

if nTR>5 && hs.countDown.Running == "off" % work for EPI without tricks
    hs.countDown.UserData = datetime(nam.date) + seconds((nTR-1)*asc_header(s,'alTR[0]')/1e6);
    start(hs.countDown);
end
try stopAt = str2double(regexp(s.ImageComments, '(?<=stopAt:)\d*', 'match', 'once')); 
catch, stopAt = inf;
end

try thk = s.SpacingBetweenSlices; catch, thk = s.SliceThickness; end
img = double(img(:,:,:));
img0 = img;
p = refVol(img, [s.PixelSpacing' thk]);
hs.img.UserData(1) = mean(img(:));
viewer = findall(0, 'Type', 'figure', 'Tag', 'nii_viewer');
if nTR<6 && isempty(viewer), overlay(hs.fig); end

R1 = inv(p.R0);
m6 = zeros(2,6);
dict = dicm_dict('', {'Rows' 'Columns' 'BitsAllocated' 'InstanceNumber'});
for i = 2:nTR
    tEnd = datetime + seconds(20); % treat as stopped series if timeout
    iNam = sprintf('%s%06i_*.dcm', bNam, nTE*nCh*(i-1)+1);
    while datetime<tEnd
        nam = dir(iNam);
        if isempty(nam), pause(0.2); else, break; end
        if startsWith(hs.MMSS.String{2}, 'Finished'), return; end
        if hs.countDown.Running=="on", countDown(hs.countDown, 0, hs); end
        if isfield(hs, 'serial'), serialRead(hs.serial); end
    end
    if isempty(nam), return; end % timeout
    s = dicm_hdr_wait(nam, dict);
    img = dicm_img(s);
    if ~isequal(numel(img), numel(img0)), return; end % bad last Instance
    img = double(img(:,:,:));
    hs.img.UserData(i) = mean(img(:)); 
    set_img(hs.img, img); hs.instnc.String = num2str(i);
    hs.slider.Value = i; % show progress
    if isDTI, continue; end % give up for now
    

    p.F.Values = smooth_mc(img, p.sz);
    [m6(2,:), R1] = moco_estim(p, R1);
    a = abs(m6(2,:) - m6(1,:)); m6(1,:) = m6(2,:);
    hs.fd.YData(i) = sum([a(1:3) a(4:6)*50]); % 50mm: brain radius
    
    a = img(:) - img0(:);
    hs.dv.YData(i) = sqrt(a'*a / numel(a)) / p.mean;
    img0 = img;
    if isfield(hs, 'serial') && hs.serial.UserData.send
        hs.serial.write(hs.dv.YData(i)/hs.ax.YAxis(1).Limits(2)*255, 'uint8');
    end
    
    a = hs.dv.YData(1:i); a = [0 a(~isnan(a))];
    dy = hs.ax.UserData; fd = hs.fd.YData(1:i);
    N = {numel(a) sum(a<dy(2)) sum(a<dy(3))};
    hs.table.Data(1,3:6) = [N mean(fd, 'omitnan')];
    for j = 1:2, hs.pct(j).String = sprintf('%.3g%%', N{j+1}/N{1}*100); end
    if i==2, set(hs.pct, 'Visible', 'on'); end
    if N{2}>=stopAt, [~, ~] = system(['touch "' hs.rootDir 'StopScan"']); end
    % drawnow; % update instance for offline test
end

%% new series or new subj: update hs.subj and hs.series.UserData
% XA30 dicm name format: 002_000003_000001_1..UID.dcm
% 1st for Study, 2nd for Series, 3rd for Instance, then long UID, .dcm
% Only Instance starts with 1 and continuous (the code relies on this)
function new = new_series(hs)
try setCountDown(hs); end
f = hs.subj.UserData;
if ~isempty(f) % check next series for current subj
    nams = dir([f '/*_*_000001*.dcm']); % 1st instance for all series
    [~, a] = sort([nams.datenum]); nams = nams(a);
    i = find(startsWith({nams.name}, hs.series.UserData), 1);
    if i < numel(nams)
        hs.series.UserData = nams(i+1).name(1:11);
        new = true; return;
    end
end
new = false;
fs = dir([hs.rootDir '20*']); % subj folder format: yyyymmdd.PatientName.PatientID
valid = cellfun(@(c)~isempty(regexp(c,'\d{8}\.[\w\d]+\.[\w\d]+','once')), {fs.name});
fs = fs(valid & [fs.isdir]);
for i = numel(fs):-1:1
    subj = regexp(fs(i).name, '(?<=\d{8}\.).*?(?=\.)', 'match', 'once');
    subj = regexprep(subj, '[_\s]', '');
    if isfile([hs.logDir subj '.mat']), continue; end
    closeSubj(hs.fig);
    hs.subj.UserData = fullfile(fs(i).folder, fs(i).name);
    hs.subj.String = subj;
    nams = dir([hs.subj.UserData '/*_*_000001*.dcm']);
    [~, a] = sort([nams.datenum]); nams = nams(a);
    hs.series.UserData = nams(1).name(1:11);
    new = true; return;
end
if ~isfile([hs.rootDir 'dClock']); return; end % rest only for RTMM computer
QC_report(hs.subj);

% Backup subj folder around 3AM
if ~isfolder(hs.backupDir) || abs(datetime-datetime('today')-hours(3))>seconds(9); return; end
fs = fs(datetime-datetime({fs.date})>days(2));
for i = 1:numel(fs)
    if isfolder([hs.backupDir fs(i).name]), continue; end
    movefile([hs.rootDir fs(i).name], hs.backupDir);
end

%% Initialize GUI for a new series
function init_series(hs, s, nTR)
fid = fopen([hs.rootDir 'currentSeries.txt'], 'w');
fprintf(fid, '%s_%s_%s', s.PatientName, asc_header(s, 'tProtocolName'), s.AcquisitionDateTime(3:12));
fclose(fid);

set(hs.slider, 'Max', nTR, 'Value', 1, 'UserData', seriesBase(s.Filename));
if nTR==1, hs.slider.Visible = 'off';
else, set(hs.slider, 'SliderStep', [1 1]./(nTR-1)); hs.slider.Visible = 'on';
end
hs.dv.YData = nan(nTR,1); hs.fd.YData = nan(nTR,1); hs.img.UserData = nan(nTR,1);
if nTR>1
    set([hs.ax hs.pct], 'Visible', 'on');
    hs.ax.XLim(2) = nTR + 0.5;
else
    set([hs.ax hs.pct], 'Visible', 'off');
end
hs.resp(1).Parent.XLim(2) = hs.ax.XLim(2);
set(hs.resp, 'XData', nan, 'YData', 1); update_resp(hs);
set([hs.instnc hs.pct], 'String', '');
figure(hs.fig); drawnow; % bring GUI front if needed
if contains(s.ImageType, '\MOCO'), return; end
hs.table.Data = [{s.SeriesDescription s.SeriesNumber 1 [] [] []}; hs.table.Data];
hs.fig.UserData.hdr{end+1} = s; % 1st instance

pat = asc_header(s,'sPat.lAccelFactPE', 1);
if pat==1, thr = 0.06; else, thr = 0.08; end % arbitrary
h = findobj(hs.fig, 'Type', 'uimenu', 'Label', '&DVARS Threshold');
thrs = str2double(get(h.Children, 'Label'));
[~, i] = min(abs(thrs-thr));
DV_yLim(h.Children(i));

%% return series base name from a file name
function nam = seriesBase(fname)
[pth, nam] = fileparts(fname);
nam = fullfile(pth, nam(1:11));

%% Set img and img axis
function set_img(h, img)
img = img(:,:,:);
d = size(img);
if numel(d)>2
    sz = getpixelposition(h.Parent);
    nMos = ceil(max(sz(3:4)*2./d(1:2))).^2;
    img(:,1,:) = 0; img(1,:,:) = 0; % visual divider
    img = vol2mos(img(:,:,unique(round(linspace(1, d(3), nMos)))));
    d = size(img);
end
set(h.Parent, 'CLim', [0 imgClim(img)], 'XLim', [0 d(2)]+0.5, 'YLim', [0 d(1)]+0.5);
h.CData = img;

%% get some series information
function c = seriesInfo(s)
c{1} = s.SeriesDescription;
if numel(c{1})>24, c{1} = [c{1}(1:16) '...' c{1}(end-3:end)]; end
c{2} = sprintf('Series %g', s.SeriesNumber);
% if s.StudyID~="1", c{2} = ['Study ' s.StudyID ', ' c{2}]; end
c{3} = dicmDT(s, 'h:mm:ss a');
c{4} = dicmDT(s, 'eee, MMM d, y');
try c{5} = sprintf('TR = %g', asc_header(s, 'alTR[0]')/1000); end

%% toggle FD display on/off
function toggleFD(h, ~)
hs = guidata(h);
if hs.fd.Visible == "on"
    set([hs.fd hs.ax.YAxis(2)], 'Visible', 'off');
    h.Label = 'Show FD plot';
else
    set([hs.fd hs.ax.YAxis(2)], 'Visible', 'on');
    h.Label = 'Hide FD plot';
end

%% Set FD plot y-axis limit
function FD_range(h, ~)
hs = guidata(h);
dy = str2double(h.Label) * (0:3) / 3;
yyaxis(hs.ax, 'right'); set(hs.ax, 'YTick', dy, 'YLim', dy([1 4]));

%% Set DVARS plot y-axis limit, and update table
function DV_yLim(h, ~)
hs = guidata(h);
dy = str2double(h.Label) * (0:3);
hs.ax.UserData = dy;
for i = 1:numel(hs.fig.UserData.DV)
    a = hs.fig.UserData.DV{i}; a = [0 a(~isnan(a))];
    if numel(a)<2, continue; end
    hs.table.Data(end-i+1, 4:5) = {sum(a<dy(2)) sum(a<dy(3))};
end
yyaxis(hs.ax, 'left'); set(hs.ax, 'YTick', dy, 'YLim', dy([1 4]));
rect = findobj(hs.ax, 'type', 'Rectangle');
for i = 1:3, rect(i).Position([2 4]) = dy([i 2]); end
DV = hs.dv.YData; a = [0 DV(~isnan(DV))];
for i = 1:2, hs.pct(i).String = sprintf('%.3g%%', sum(a<dy(i+1))/numel(a)*100); end

%% Table-click callback: show moco/series info and image if avail
function tableCB(h, evt)
if isempty(evt.Indices) || evt.Indices(1,2)>2, return; end
hs = guidata(h);
C = h.Data;
iT = evt.Indices(1,1);
iR = size(C,1) - iT + 1;
dV = hs.fig.UserData.DV{iR};
hs.fd.YData = hs.fig.UserData.FD{iR};
hs.dv.YData = dV;
nTR = numel(dV);
if numel(dV)<2 || isnan(dV(2))
    set([hs.ax hs.pct], 'Visible', 'off');
else    
    set([hs.ax hs.pct], 'Visible', 'on');
    for j = 1:2, hs.pct(j).String = sprintf('%.3g%%', C{iT,j+3}/nTR*100); end
    hs.resp(1).Parent.XLim(2) = hs.ax.XLim(2);
    try
        for i = 1:3
            a = hs.fig.UserData.resp{iR}{i};
            set(hs.resp(i), 'XData', a, 'YData', ones(size(a)));
        end
    catch, set(hs.resp, 'XData', nan, 'YData', 1);
    end
    update_resp(hs);
end
hs.instnc.String = '';
hs.series.String = C{iT,1}; % in case hdr not saved
try s = hs.fig.UserData.hdr{iR}; catch, set_img(hs.img, inf(2)); return; end
if ~isfile(s.Filename), s.Filename = strrep(s.Filename, hs.rootDir, hs.backupDir); end

iIN = ceil(nTR/2); % start with middle Instance if avail
nam = dir(sprintf('%s%06g.dcm', seriesBase(s.Filename), iIN));
if isempty(nam), iIN = 1; nam = s;
else, nam = fullfile(nam(1).folder, nam(1).name);
end
hs.instnc.String = num2str(iIN);
set(hs.slider, 'Max', nTR, 'Value', iIN, 'UserData', seriesBase(s.Filename));
if nTR == 1, hs.slider.Visible = 'off';
else, set(hs.slider, 'SliderStep', [1 1]./(nTR-1)); hs.slider.Visible = 'on';
end
hs.ax.XLim(2) = numel(hs.dv.YData) + 0.5;
hs.series.String = seriesInfo(s);
try img = dicm_img(nam); catch, img = ones(2)*0.94; end
set_img(hs.img, img);

%% Load subj data to review
function loadSubj(h, ~)
hs = guidata(h);
[fname, pName] = uigetfile([hs.logDir '*.mat'], 'Select MAT file for a Patient');
if isnumeric(fname), return; end
load([pName '/' fname], 'T3'); % T4 later
hs.fig.UserData = T3.Properties.UserData; 
DV = hs.fig.UserData.DV;
[N, M] = size(T3);
C = flip(table2cell(T3), 1); C(:,6) = C(:,M);
dy = hs.ax.UserData;
for i = 1:N
    a = DV{N-i+1}; a = [0 a(~isnan(a))];
    if numel(a)<2, C(i,4:5) = {[] []};
    else, C(i,4:5) = {sum(a<dy(2)) sum(a<dy(3))};
    end
    if M<4, C{i,3} = numel(a); end % did not save this before 240205
end 
hs.table.Data = C;
fnam = hs.fig.UserData.hdr{end}.Filename;
[hs.subj.UserData, nam] = fileparts(fnam);
subj = regexp(hs.subj.UserData, '(?<=\d{8}\.).*?(?=\.)', 'match', 'once');
hs.subj.String = regexprep(subj, '[_\s]', '');
hs.series.UserData = nam(1:11);
tableCB(hs.table, struct('Indices', [1 1])); % show top series

%% close subj
function closeSubj(h, ~)
hs = guidata(h);
hs.table.Data = cell(0,6);
hs.img.CData = ones(2)*0.94;
hs.fd.YData = nan; hs.dv.YData = nan;
c = {{}};
hs.fig.UserData = struct('FD', c, 'DV', c, 'hdr', c, 'resp', c, 'GM', c);
set([hs.subj hs.series hs.instnc], 'String', '');
set([hs.subj hs.series], 'UserData', '');
set(hs.pct, 'Visible', 'off');
close(findall(0, 'Type', 'figure', 'Tag', 'nii_viewer')); % last subj if any

%% Re-do current subj: useful in case of error during a session
function redoSubj(h, ~)
hs = guidata(h);
subj = hs.subj.String;
if isempty(subj), return; end
if ~isfolder(hs.subj.UserData)
    fprintf(2, 'Image for %s deleted?\n', subj);
    return;
end
try delete([hs.logDir subj '*.mat']); end
closeSubj(hs.fig)
hs.subj.String = subj;

%% Get reference vol info. Adapted from nii_moco.m
function p = refVol(img, pixdim)
d = size(img);
p.R0 = diag([pixdim 1]); % no need for real xform_mat here
p.R0(1:3, 4) = -pixdim .* (d/2); % make center voxel [0 0 0]

sz = pixdim;
if all(abs(diff(sz)/sz(1))<0.05) && sz(1)>2 && sz(1)<4 % 6~12mm
    sz = 3; % iso-voxel, 2~4mm res, simple fast smooth
else
    sz = 9 ./ sz'; % 9 mm seems good
end

% resample ref vol to isovoxel (often lower-res)
d0 = d-1;
dd = 4 ./ pixdim; % use 4 mm grid for alignmen
[i, j, k] = ndgrid(0:dd(1):d0(1)-0.5, 0:dd(2):d0(2)-0.5, 0:dd(3):d0(3)-0.5);
I = [i(:) j(:) k(:)]';
a = rng('default'); I = I + rand(size(I))*0.5; rng(a); % used by spm
V = smooth_mc(img, sz);
F = griddedInterpolant({0:d0(1), 0:d0(2), 0:d0(3)}, V, 'linear', 'none');
V0 = F(I(1,:), I(2,:), I(3,:)); % ref: 1 by nVox
I(4,:) = 1; % 0-based ijk: 4 by nVox
I = p.R0 * I; % xyz of ref voxels

% compute derivative to each motion parameter in ref vol
dG = zeros(6, numel(V0));
dd = 1e-6; % delta of motion parameter, value won't affect dG much
R0i = inv(p.R0); % speed up a little
for i = 1:6
    p6 = zeros(6,1); p6(i) = dd; % change only 1 of 6
    J = R0i * rigid_mat(p6) * I; %#ok<*MINV>
    dG(i,:) = F(J(1,:), J(2,:), J(3,:)) - V0; % diff now
end
dG = dG / dd; % derivative

% choose voxels with larger derivative for alignment: much faster
a = sum(dG.^2); % 6 derivatives has similar range
ind = a > std(a(~isnan(a)))/10; % arbituray threshold. Also exclude NaN
p.dG = dG(:, ind);
p.V0 = V0(ind);
p.mean = mean(p.V0);
p.mm = I(:, ind);
F.GridVectors = {0:d(1)-1, 0:d(2)-1, 0:d(3)-1};
p.F = F;
p.sz = sz;

%% motion correction to ref-vol. From nii_moco.m
function [m6, R] = moco_estim(p, R)
if isnan(R(1)), R = inv(p.R0); end
for iter = 1:64
    J = R * p.mm; % R_rst*J -> R0*ijk
    V = p.F(J(1,:), J(2,:), J(3,:));
    ind = ~isnan(V); % NaN means out of range
    dV = p.V0(ind) - V(ind);
    if sum(ind)<32 || dV*dV'/numel(dV)>p.mean*p.mean, R = nan(4); break; end
    a = p.dG(:, ind);
    b = (a * a') \ (a * dV'); % dG(:,ind)'\dV' estimate p6 from current R
    R = R * rigid_mat(b); % inv(inv(rigid_mat(p6)) * inv(R_rst))
    if b'*b < 1e-4, break; end % little effect, stop
end

M = p.R0 * R; % inv(R_rst / Rref)
m6 = -[M(1:3, 4)' atan2(M(2,3), M(3,3)) asin(M(1,3)) atan2(M(1,2), M(1,1))];

%% Translation (mm) and rotation (deg) to 4x4 R. Order: ZYXT
function R = rigid_mat(p6)
ca = cosd(p6(4:6)); sa = sind(p6(4:6));
rx = [1 0 0; 0 ca(1) -sa(1); 0 sa(1) ca(1)]; % 3D rotation
ry = [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)];
rz = [ca(3) -sa(3) 0; sa(3) ca(3) 0; 0 0 1];
R = rx * ry * rz;
R = [R p6(1:3); 0 0 0 1];

%% Simple gaussian smooth for motion correction, sz in unit of voxels
function out = smooth_mc(in, sz)
out = double(in);
if numel(unique(in))<5, return; end
if all(abs(diff(sz)/sz(1))<0.05) && abs(sz(1)-round(sz(1)))<0.05 ...
        && mod(round(sz(1)),2)==1
    out = smooth3(out, 'gaussian', round(sz)); % sz odd integer
    return; % save time for special case
end

d = size(in);
I = {1:d(1) 1:d(2) 1:d(3)};
n = sz/3;
if numel(n)==1, n = n*[1 1 1]; end
J = {1:n(1):d(1) 1:n(2):d(2) 1:n(3):d(3)};
intp = 'linear';
F = griddedInterpolant(I, out, intp);
out = smooth3(F(J), 'gaussian'); % sz=3
F = griddedInterpolant(J, out, intp);
out = F(I);

%% Subfunction: get a parameter in CSA series ASC header: MrPhoenixProtocol
function val = asc_header(s, key, dft)
if nargin>2, val = dft; else, val = []; end
csa = 'CSASeriesHeaderInfo';
if ~isfield(s, csa) % in case of multiframe
    try s.(csa) = s.SharedFunctionalGroupsSequence.Item_1.(csa).Item_1; end
end
if isfield(s, 'Private_0029_1020') && isa(s.Private_0029_1020, 'uint8')
    str = char(s.Private_0029_1020(:)');
    str = regexp(str, 'ASCCONV BEGIN(.*)ASCCONV END', 'tokens', 'once');
    if isempty(str), return; end
    str = str{1};
elseif isfield(s, 'MrPhoenixProtocol') % X20A
    str = s.MrPhoenixProtocol;
elseif ~isfield(s, csa), return; % non-siemens
elseif isfield(s.(csa), 'MrPhoenixProtocol') % most Siemens dicom
    str = s.(csa).MrPhoenixProtocol;
elseif isfield(s.(csa), 'MrProtocol') % older version dicom
    str = s.(csa).MrProtocol;
else, return;
end

% tSequenceFileName  = ""%SiemensSeq%\gre_field_mapping""
expr = ['\n' regexptranslate('escape', key) '\s*=\s*(.*?)\n'];
str = regexp(str, expr, 'tokens', 'once');
if isempty(str), return; end
str = strtrim(str{1});

if strncmp(str, '""', 2) % str parameter
    val = str(3:end-2);
elseif strncmp(str, '"', 1) % str parameter for version like 2004A
    val = str(2:end-1);
elseif strncmp(str, '0x', 2) % hex parameter, convert to decimal
    val = sscanf(str(3:end), '%x', 1);
else % decimal
    val = sscanf(str, '%g', 1);
end

%% Wait till a file copy is complete
function s = dicm_hdr_wait(nam, dict)
tEnd = datetime + seconds(2); % wait till file ready
nam = [nam.folder '/' nam.name];
while 1
    s = dicm_hdr(nam, dict);
    try
        if s.PixelData.Start+s.PixelData.Bytes<=s.FileSize || datetime>tEnd
            return;
        end
    catch me
        if datetime>tEnd, fprintf(2,'%s\n', nam); rethrow(me); end % give up
    end
    pause(0.1);
end

%% User closing GUI: stop and delete timer
function closeFig(fh, ~)
hs = guidata(fh);
delete(fh);
try hs.timer.StopFcn = ''; end
try tObj = timerfindall; stop(tObj); delete(tObj); end
clear hs; % close serialport

%% menu callback for both DERIVED and _SBRef
function toggleChecked(h, ~)
if h.Checked == "on", h.Checked = 'off'; else, h.Checked = 'on'; end

%% Increase/Decrease image CLim
function setCLim(h, ~)
hs = guidata(h);
ax = hs.img.Parent;
if startsWith(h.Label, 'Increase'), ax.CLim(2) = ax.CLim(2)*0.8;
elseif startsWith(h.Label, 'Decrease'), ax.CLim(2) = ax.CLim(2)*1.2;
else, ax.CLim(2) = imgClim(hs.img.CData);
end

%% show series in nii_viewer
function view_3D(h, ~)
hs = guidata(h);
nams = dir([hs.slider.UserData '*.dcm']);
if isempty(nams), return; end
nams = strcat(nams(1).folder, '/', {nams.name});
nii = dicm2nii(nams, ' ', 'no_save');
nii_viewer(nii);

%% overlay series onto T1w or Scout if avail
function overlay(h, ~)
hs = guidata(h);
hdrs = hs.fig.UserData.hdr;
if isempty(hdrs), return; end
is3D = @(c)c.MRAcquisitionType=="3D" && ~contains(c.ImageType, 'DERIVED');
is3D = cellfun(is3D, hdrs);
if ~any(is3D), is3D = cellfun(@(c)c.MRAcquisitionType=="3D", hdrs); end
if sum(is3D)>1
    isT1 = cellfun(@(c)contains(c.SequenceName, 'fl3d1'), hdrs);
    isT1 = isT1 & is3D;
    if any(isT1), is3D = isT1; end
end
if ~any(is3D), view_3D(h); return; end % no T1, just show in nii_viewer
is3D = find(is3D, 1, 'last');
a = seriesBase(hdrs{is3D}.Filename);
nams = dir([a '*.dcm']);
if isempty(nams), nams = dir([strrep(a, hs.rootDir, hs.backupDir) '*.dcm']); end
nams = strcat(nams(1).folder, '/', {nams.name});
T1w = dicm2nii(nams, ' ', 'no_save');
nams = dir([hs.slider.UserData '*.dcm']);
nams = strcat(nams(1).folder, '/', {nams.name});
epi = dicm2nii(nams, ' ', 'no_save');
fh = nii_viewer(T1w, epi); fh.Position(2) = fh.Position(2)-200;
nii_viewer('LocalFunc', 'nii_viewer_cb', [], [], 'center', fh);

%% slider callback: show img if avail
function sliderCB(h, ~)
hs = guidata(h);
if isempty(h.UserData), return; end
h.Value = round(h.Value);
base = sprintf('%s%06u*.dcm', h.UserData, h.Value);
nam = dir(base);
if isempty(nam), nam = dir(strrep(base, hs.rootDir, hs.backupDir)); end
if isempty(nam), return; end
set_img(hs.img, dicm_img(fullfile(nam(1).folder, nam(1).name)));
hs.instnc.String = num2str(h.Value);

%% Timer StopFunc: Save result, start timer with delay, even after error
function saveResult(obj, ~)
hs = guidata(obj.UserData);
set([hs.menu hs.table hs.slider], 'Enable', 'on');
toSave = false;
if size(hs.table.Data,1)
    N = numel(dir([seriesBase(hs.fig.UserData.hdr{end}.Filename) '*.dcm']));
    if N > hs.table.Data{1,3}
        hs.table.Data{1,3} = N; % update instances
        toSave = true;
    end
end
if size(hs.table.Data,1) > numel(hs.fig.UserData.FD) % new series done
    if isfield(hs, 'serial'), hs.serial.UserData.send = false; end % stop sending
    hs.fig.UserData.FD{end+1} = hs.fd.YData;
    hs.fig.UserData.DV{end+1} = hs.dv.YData;
    hs.fig.UserData.GM{end+1} = hs.img.UserData;
    hs.fig.UserData.resp{end+1} = {hs.resp.XData};
    toSave = true;
end
if toSave
    T3 = cell2table(flip(hs.table.Data(:,[1:3 6]), 1), ...
        'VariableNames', {'Description' 'SeriesNumber' 'Instances' 'MeanFD'});
    T3.Properties.UserData = hs.fig.UserData;
    save([hs.logDir hs.subj.String], 'T3');
end
if new_series(hs), obj.StartDelay = 0.1; else, obj.StartDelay = 5; end
start(obj);

%% Serial BytesAvail callback: update response: 1=missed, 2=incorrect, 3=correct 
function serialRead(s, ~)
if s.NumBytesAvailable<1, return; end
b = s.read(1, 'uint8');
hs = guidata(s.UserData.fig);
if     b == '?', s.write('RTMM', 'char'); % identity
elseif b == 'P', s.write(hs.subj.String, 'char'); % PatientName
elseif b == 'T', s.write(strjoin(hs.MMSS.String), 'char'); % count down
elseif b == 'M', s.UserData.send = true; % start to send motion info
elseif b == 'Q' % stim computer asks to stop scan
    [~, ~] = system(['touch "' hs.rootDir 'StopScan"']);
elseif b>0 && b<4 % response
    x = find(~isnan(hs.dv.YData), 1, 'last');
    if x == numel(hs.dv.YData), return; end % not running
    if isempty(x), x = 0; end
    hs.resp(b).XData(end+1) = x + 1; hs.resp(b).YData(end+1) = 1;
    update_resp(hs);
end

%% update response text
function update_resp(hs)
n = cellfun(@numel, {hs.resp.XData}) - 1;
h = hs.resp(1).Parent.Title;
if ~any(n>0), h.Visible = 'off'; return; end
h.Visible = 'on';
h.String = "Missed " + n(1) + ", \color{red}Incorrect " + n(2) + ...
    ", \color[rgb]{0 0.8 0}Correct " + n(3) + ", \color{blue}Total " + sum(n);

%% start countdown
function setCountDown(hs)
nam = [hs.rootDir 'SyngoMeas'];
if ~isfile(nam), return; end
c0 = fileread(nam); pause(0.2); c = fileread(nam);
if ~isequal(c0, c), pause(1); c = fileread(nam); end
delete(nam); 
% From scanner: "RunStartTime" "ProtocolName" TotalScanTimeSec
c = regexp(c, '"(.*?)" "(.*?)" (\d+)', 'tokens', 'once');
if isempty(c), return; end % MeasFinished?
tStart = datetime(c{1}, 'Format', 'yyyy/MM/dd-HH:mm:ss.SSSSSS');
dClock = load([hs.rootDir 'dClock'], '-ascii');
tFnsh = tStart + seconds(str2double(c{3})-dClock);
if tFnsh-datetime < seconds(2), return; end
hs.countDown.UserData = tFnsh;
if numel(c{2})>24, c{2} = [c{2}(1:16) '...' c{2}(end-3:end)]; end
hs.MMSS.String = {c{2} ''};
if hs.countDown.Running=="off", start(hs.countDown); end

%% timer func to show scanning time
function countDown(tObj, ~, hs)
t = tObj.UserData - datetime; t.Format = 'mm:ss';
if t<seconds(1), stop(tObj); hs.MMSS.String = {'' ''}; return; end
hs.MMSS.String{2} = ['Scanning ' char(t)];
nam = hs.rootDir + "SyngoMeas";
if ~isfile(nam), return; end
if ~startsWith(fileread(nam), "Finished"), return; end
delete(nam);
stop(tObj);
hs.MMSS.String{2} = ['Finished ' char(t)];

%% get CLim for dicom img
function mx = imgClim(img)
im = double(img(:));
im = im(im>max(im)/10);
mx = mean(im) + 2*std(im);

%% Create QC report
function QC_report(h, ~)
hs = guidata(h);
if nargin>1 % calling from menu
    subj = hs.subj.String;
else
    nam = dir([hs.rootDir 'closed_*']);
    if isempty(nam) || datetime-datetime(nam(1).date)<seconds(1) || ...
            isfile([hs.rootDir 'EyelinkRecording.mat']); return; end
    nam = [hs.rootDir nam(1).name];
    done = onCleanup(@()movefile(nam, strrep(nam, 'closed_', 'done_')));
    subj = regexp(nam, '(?<=closed_)\d{4,6}\w{2}$', 'match', 'once');
end
rmQC = onCleanup(@()delete('./tmp_QC_*.pdf'));
try load([hs.rootDir 'RTMM_log/' subj '.mat'], 'T3'); catch, return; end
uDat = T3.Properties.UserData;

delete('./tmp_QC_*.pdf');
fig = figure('Position', [10 30 [8.5 11]*96], 'Units', 'normalized');
closeFig = onCleanup(@()delete(fig));
set(fig, 'Color', 'w', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
ax = axes(fig, 'Position', [0.01 0.955 0.98 0.045], 'Visible', 'off');
try imshow('./logo.png', 'Parent', ax); ax.HandleVisibility = 'off'; end

layout = getpref('nii_viewer_para', 'layout');
if layout ~= 1
    setpref('nii_viewer_para', 'layout', 1);
    cln = onCleanup(@()setpref('nii_viewer_para', 'layout', layout));
end

ax = axes(fig, 'Position', [0.1 0.92 0.8 0.03], 'Visible', 'off');
text(ax, 0.5, 1, subj, 'FontSize', 18, 'HorizontalAlignment', 'center');
s = uDat.hdr{1};
text(ax, 0.5, 0, dicmDT(s,'eeee MMM d, y'), 'FontSize', 12, 'HorizontalAlignment', 'center');
f = [fileparts(s.Filename) filesep];
nam1 = dir([f '*_*_000001*.dcm']); % 1st instance for all series
[~, a] = sort([nam1.datenum]); nam1 = {nam1(a).name};
tbl = cell(numel(nam1), 5);
dict = dicm_dict('', {'AcquisitionDateTime' 'SeriesNumber' 'SeriesDescription'});
for i = 1:size(tbl,1)
    nams = dir([f nam1{i}(1:11) '*.dcm']);
    s = dicm_hdr([f nams(1).name], dict);
    ind = T3.SeriesNumber==s.SeriesNumber & strcmp(s.SeriesDescription, T3.Description);
    if any(ind), a = T3.MeanFD{find(ind,1)}; else, a = []; end
    tbl(i,:) = {s.SeriesNumber dicmDT(s,'HH:mm:ss') numel(nams) s.SeriesDescription a};
end
tbl = cellfun(@num2str, tbl, 'UniformOutput', false); % for left-align
vName = {'SeriesNumber' 'Time' 'TotalInstances' 'Description' 'meanFD'};
h = uitable(fig, 'Units', 'normalized', 'Position', [0.04 0.04 0.9 0.85], ...
    'FontSize', 12, 'RowName', [], 'ColumnName', vName, 'Data', tbl, ...
    'ColumnWidth', {96 96 108 342 82});
i = size(tbl, 1);
if i>40, h.FontSize = max(6, fix(450/i)); end
warning('off', 'MATLAB:print:ExcludesUIInFutureRelease');
newPage(fig); y = 0.95;

el = {}; clear EL;
nam = dir([hs.rootDir 'RTMM_log/' subj '_*.edf']);
while ~isempty(nam)
    tm = []; tg = []; pa = [];
    for i = numel(nam):-1:1 % in case of multiple files
        [~, a] = evalc("edfmex('"+nam(i).folder+"/"+nam(i).name+"')");
        tm = [a.FSAMPLE.time tm]; %#ok
        tg = [a.FSAMPLE.buttons tg]; %#ok
        pa = [a.FSAMPLE.pa pa]; %#ok
    end
    if numel(unique(pa(1,:)))>9, pa = pa(1, :);
    elseif numel(unique(pa(2,:)))>9, pa = pa(2, :);
    else, break;
    end
    i = find(arrayfun(@(c)startsWith(char(c.message), 'DICM_'), a.FEVENT), 1, 'last');
    a = double(tm)/1000 - double(a.FEVENT(i).sttime)/1000 + ...
        sscanf(a.FEVENT(i).message, 'DICM_secs=%g');
    EL.Hz = 100; % resample pa to 100 Hz
    EL.t0 = a(1); % recording start time in Syngo secs of the day
    EL.tg = a(diff([0 bitget(tg, 5)])>0); % trigger time
    EL.pa = interp1(a, pa/55^2, a(1):1/EL.Hz:a(end)); % convert roughly to degree
    break;
end

for i = 1:numel(uDat.hdr)
    s = uDat.hdr{i};
    series = sprintf('%s (Series %g)', s.SeriesDescription, s.SeriesNumber);
    if contains(s.SequenceName, {'epfid2d' 'mbPCASL' 'ep_b0'}) % EPI/ASL/Diff
        fd = uDat.FD{i};
        N = find(~isnan(fd), 1, 'last');
        if isempty(N) || N<11 || all(fd==0), continue; end % skip slice check etc
        if y<0.35, newPage(fig); y = 0.95; end
        if exist('EL', 'var') && contains(s.SequenceName, 'epfid2d')
            t0 = seconds(timeofday(dicmDT(s)));
            [mi, j] = min(abs(EL.tg - t0));
            if mi<5 % recording started late if mi too large
                TR = asc_header(s, 'alTR[0]')*1e-6;
                while j>1 && abs(diff(EL.tg([j-1 j]))-TR)<0.1, j = j-1; end
                % fprintf('%4.2f\n', EL.tg(j)-t0); % normally <1
                j = floor((EL.tg(j)-EL.t0) * EL.Hz);
                nP = round(N * TR * EL.Hz);
                el{end+1} = {EL.pa(j+(1:nP)) series}; %#ok save for users
                y = y - 0.4;
                ax0 = axes(fig, 'Position', [0.1 y+0.27 0.8 0.08]);
                plot(ax0, el{end}{1}, 'b');
                set(ax0, 'XLim', [1 nP], 'XTick', []);
                ylabel(ax0, 'Pupil Size', 'Color', 'b');
            end
        else, y = y - 0.32; clear ax0;
        end
        ax = axes(fig, 'Position', [0.1 y+0.19 0.8 0.08]);
        gm = uDat.GM{i}(1:N);
        gm = (gm/mean(gm) - 1)*100; rg = ceil(std(gm)); % global mean
        % gm = [0; diff(gm)]/mean(gm) * 100; % delta GM
        plot(ax, gm, '.-m');
        set(ax, 'YTick', [], 'YLim', [-rg rg], 'XLim', [0 N+1], 'XTick', []);
        ylabel(ax, 'GM (%)', 'Color', 'm');
        try ax = ax0; end
        title(ax, series, 'Interpreter', 'none', 'FontSize', 12);
        % ylabel(ax, [char(916) 'GM (%)']);
        ax = axes(fig, 'Position', [0.1 y+0.03 0.8 0.16]);
        yyaxis(ax, 'right'); plot(ax, fd(1:N), '.-');
        ylabel(ax, 'FD (mm)'); xlabel(ax, 'Volume Number');
        set(ax, 'YTick', 0:0.6:1.8, 'YLim', [0 1.8], 'XLim', [0 N+1]);
        str = sprintf('meanFD=%.2g', mean(fd(2:N)));
        text(ax, 0.8, 0.9, str, 'Units', 'normalized', 'Color', [0.85 0.32 0.1]);
        yyaxis(ax, 'left'); plot(ax, uDat.DV{i}(1:N), '.-');
        ylabel(ax, 'DVARS'); set(ax, 'YTick', 0:0.12:0.36, 'YLim', [0 0.3]);
    elseif contains(s.SequenceName, {'tfl3d' 'tfl_me3d' 'spc' 'tse2d' 'fm2d' 'epse2d'})
        if endsWith(s.SeriesDescription, '_setter'), continue; end
        nams = dir([seriesBase(s.Filename) '000001*.dcm']);
        nams = strcat(nams(1).folder, '/', {nams.name});
        nii = dicm2nii(nams, ' ', 'no_save');
        if isempty(nii), continue; end
        fh = nii_viewer(nii);
        nii_viewer('LocalFunc', 'nii_viewer_cb', [], [], 'center', fh);
        hsV = guidata(fh);
        drawnow; F = getframe(fh, hsV.frame.Position); close(fh); img = F.cdata;
        sz = size(img); sz = sz([2 1]) ./ [8.5 11]/96;
        if sz(1)>0.8, sz = sz / sz(1) * 0.8; end
        y = y - sz(2) - 0.08;
        if y<0.02 && y+sz(2)*0.2>0.02, sz = sz * 0.8; y = y + 0.2*sz(2); end
        if y<0.02, newPage(fig); y = 0.92-sz(2); end
        ax = axes(fig, 'Position', [(1-sz(1))/2 y+0.02 sz], 'Visible', 'off');
        imshow(img, 'Parent', ax);
        title(ax, series, 'Interpreter', 'none', 'FontSize', 12);
    end
end
newPage(fig); close(fig);

if ~isempty(el), save([hs.rootDir 'RTMM_log/' subj '_eye.mat'], 'el'); end
pdfNam = hs.rootDir+"RTMM_log/"+subj+"_"+dicmDT(s,'yyMMdd')+"_QC.pdf";
PDFs = dir('./tmp_QC_*.pdf');
mergePdfs({PDFs.name}, pdfNam); 

%% export fig to a new PDF, called by QC_report()
function newPage(fig)
nam = dir('./tmp_QC_*.pdf');
if isempty(nam), i = 1; else, i = str2double(nam(end).name(8:10)) + 1; end
print(fig, sprintf('./tmp_QC_%03i.pdf', i), '-dpdf');
% exportapp(fig, sprintf('./tmp_QC_%03i.pdf', i)); % for uifigure
clf(fig);

%% merge PDFs, by Benjamin Großmann
% https://www.mathworks.com/matlabcentral/fileexchange/89127-merge-pdf-documents
function mergePdfs(INs, out)
memSet = org.apache.pdfbox.io.MemoryUsageSetting.setupMainMemoryOnly();
merger = org.apache.pdfbox.multipdf.PDFMergerUtility;
cellfun(@(f) merger.addSource(f), INs);
merger.setDestinationFileName(out);
merger.mergeDocuments(memSet);

%% 3D to mos: for display purpose
function mos = vol2mos(vol)
[nr, nc, nSL] = size(vol);
nMos = ceil(sqrt(nSL));
mos = zeros([nr nc]*nMos, 'like', vol);
for i = 0:nMos-1
    r = i*nr + (1:nr);
    for j = 1:nMos
        iSL = i*nMos + j;
        if iSL>nSL, return; end
        c = (j-1)*nc + (1:nc);
        mos(r,c) = vol(:,:,iSL);
    end
end

%% return datetime obj, or date/time str if fmt provided
function dt = dicmDT(s, fmt)
dt = datetime(s.AcquisitionDateTime, 'InputFormat', 'yyyyMMddHHmmss.SSSSSS');
if nargin>1, dt.Format = fmt; dt = char(dt); end

%%