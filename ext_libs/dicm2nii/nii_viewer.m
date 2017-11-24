function varargout = nii_viewer(fname, overlayName)
% Basic tool to visualize NIfTI images.
% 
%  NII_VIEWER('/data/subj2/fileName.nii.gz')
%  NII_VIEWER('background.nii', 'overlay.nii')
%  NII_VIEWER('background.nii', {'overlay1.nii' 'overlay2.nii'})
% 
% If no input is provided, the viewer will ask for background NIfTI. Although
% the preferred format is NIfTI, NII_VIEWER accepts any files that can be
% converted into NIfTI by dicm2nii.
% 
% Here are some features and usage.
% 
% The basic use is to open a NIfTI file to view. When a NIfTI (background) is
% open, the display always uses the image plane close to xyz axes (voxel space)
% even for oblique acquisition. The possible confusion arises if the acquisition
% was tilted with a large angle, and then the orientation labels will be less
% accurate. The benefit is that no interpolation is needed for the background
% image, and there is no need to switch coordinate system when images with
% different systems are overlaid. The display is always in correct scale at
% three axes even with non-isotropic voxels. The displayed IJK always correspond
% to left -> right, posterior -> anterior and inferior -> superior directions,
% while the NIfTI data may not be saved in this order or along these directions.
% The I-index is increasing from left to right even when the display is flipped
% as radiological convention (right on left side).
% 
% Navigation in 4D can be by mouse click, dialing IJK and volume numbers, or
% using keys (arrow keys and [ ] for 3D, and < > for volume).
% 
% After the viewer is open, dragging and dropping a NIfTI file will open it as
% background, while dropping with Ctrl key down will add it as overlay.
% 
% By default, the viewer shows full view of the background image data. The
% zoom-in always applies to three views together, and enlarges around the
% location of the crosshair. To zoom around a different location, set the
% crosshair to the interested location, and apply zoom again either by View ->
% Zoom in, or pressing Ctrl (Cmd) and +/-. View -> Set crosshair at -> center of
% view (or pressing key C) can set the crosshair to the center of display for
% all three views.
% 
% Overlays are always mapped onto the coordinate of background image, so
% interpolation (nearest/linear/cubic/spline) is usually involved. The overlay
% makes sense only when it has the same coordinate system as the background
% image, while the resolution and dimension can be different. The viewer tries
% to match any of sform and qform between the images. If there is no match, a
% warning message will show up.
% 
% A special overlay feature "Add aligned overlay" can be used to check the
% effect of registration, or overlay an image to a different coordinate system
% without creating a transformed file. It will ask for two files. The first is
% the overlay NIfTI file, and the second is either a transformation matrix file
% which aligns the overlay to the background image, or a warp file which
% transforms the overlay into the background reference.
% 
% Here is an example to check FSL alignment. From a .feat/reg folder, Open
% "highres" as background image. "Add overlay" for "example_func". If there is
% head movement between the highres and the functional scans, the overlap will
% be off. Now "Add aligned overlay" for "example_func", and use
% example_func2highres.mat as the matrix. The two dataset should overlap well if
% the alignment matrix is accurate.
% 
% Here is another example to overlay to a different coordinate system for FSL
% result. Open .feat/reg/standard.nii.gz as background. If an image like
% .feat/stats/zstat1.nii.gz is added as an overlay, a warning message will say
% inconsistent coordinate system since the zstat image is in Scanner Anat. The
% correct way is to "Add aligned overlay" for zstat1, and either use
% .feat/reg/example_func2standard.mat for linear transformation or better use
% .feat/reg/example_func2standard_warp.nii.gz if available for alignment.
% 
% When the mouse pointer is moved onto a voxel, the voxel indices, corresponding
% x/y/z coordinates and voxel value will show on the panel. If there is an
% overlay, the overlay voxel value will also show up, unless its display is off.
% When the pointer moves onto the panel or the bottom-right quadrant, the
% information for the voxel at crosshair will be displayed. The display format
% is as following with val of the top image displayed first:
%  (i,j,k)=(x,y,z): val_1 val_2 val_3
%
% Note that although the x/y/z coordinates are shared by background and overlay
% images, IJK indices are always for background image (name shown on title bar).
% 
% The mouse-over display can be turned on/off from Help -> Preferences ...
% 
% If there is a .txt label file in the same folder as the NIfTI file, like for
% AAL, the labels will be shown instead of voxel value. The txt file should have
% a voxel value and ROI label pair per line, like
%  1 Precentral_L
%  2 Precentral_R
%  3 ...
% 
% Image display can be smoothed in 3D (default is off). The smooth is slow when
% the image dimension is large, even when the current implementation of smooth
% does not consider voxel size.
% 
% Background image and overlays are listed at the top-left of the panel. All
% parameters of the bottom row of the panel are for the selected image. This
% feature is indicated by a frame grouping these parameters. Each NIfTI file has
% its own set of parameters (display min and max value, LUT, alpha, whether to
% smooth, interpolation method, and volume number) to control its display.
% Moving the mouse onto a parameter will show the meaning of the parameter.
% 
% The lastly added overlay is on the top of display and top of the file list.
% The background and overlay order can be changed by the two small buttons next
% to the list, or from Overlay -> Move selected image ...
% 
% Each NIfTI display can be turned on/off by clicking the small checkbox at the
% left side of the file (or pressing spacebar for the selected NIfTI). This
% provides a way to turn on/off an overlay to view the overlap. Most operations
% are applied to the selected NIfTI in the list, such as Show NIfTI hdr/ext
% under Window menu, Move/Close overlay under Overlay menu, and most operations
% under File menu.
% 
% A NIfTI mask can be applied to the selected image. Ideally, the mask should be
% binary, and the image corresponding to the non-zero part of the mask will be
% displayed. If non-binary mask is detected, a threshold to binarize will be
% asked. If the effect is not satisfied with the threshold, one can apply the
% mask with a different threshold without re-loading image. The way to remove a
% mask is to Remove, then Add overlay again. In case one likes to modulate the
% selected image with another NIfTI image (multiply two images), File -> Apply
% modulation will do it. If the mask image is not within 0 and 1, the lower and
% upper bound will be asked to normalize the range to [0 1]. A practical use of
% modulation is to use dti_FA map as weight for dti_V1 RGB image.
% 
% For multi-volume data, one can change the Volume Number (the parameter at
% rightmost of the panel) to check the head motion. Click in the number dialer
% and or press < or > key, to simulate movie play. It is better to open the 4D
% data as background, since it may be slower to map it to the background image.
% 
% Popular LUT options are implemented. Custom LUT can be added by File -> Load
% custom LUT. The color coding can be shown by View -> Show colorbar. There are
% several special LUTs. The "two-sided" allows to show both positive and
% negative data in one view. For example, if the display range is 3 to 10 for a
% t-map, positive T above 3 will be coded as red-yellow, and T below -3 will be
% coded as blue-green. This means the absolute display range values are used.
% 
% One of the special LUT is "lines" in red text. This is for normalized vector
% display, e.g. for diffusion vector. The viewer will refuse the LUT selection
% if the data is not normalized vector. Under this LUT, all other parameters for
% the display are ignored. The color of the "lines" is the max color of previous
% LUT. For example, if one likes to show blue vector lines, choose LUT "blue"
% first, then change it to "lines".
% 
% In case of complex image, most LUTs will display only the magnitude of the
% data, except the following three phase LUTs, where the magnitude is used as
% mask. Here is how the 3 phase LUTs encodes phase from 0 to 360 degree:
%  phase:  red-yellow monotonically,
%  phase3: red-yellow-green-yellow-red circularly, and 
%  phase6: red-yellow-green/violet-blue-cyan, with sharp change at 180 degree. 
% These LUTs are useful to visualize activation of periodic stimulus, such as
% those by expanding/shrinking ring or rotating wedge for retinotopy. To use
% this feature, one can save an complex NIfTI storing the Fourier component at
% the stimulus frequency.
% 
% Note that, for RGB NIfTI, the viewer always displays the data as color
% regardless of the LUT option. The display min and max value also have no
% effect on RGB image. There is a special LUT, RGB, which is designed to display
% non-RGB NIfTI data as RGB, e.g. FSL-style 3-volome data. 
% 
% The viewer figure can be copied into clipboard (not available for Linux) or
% saved as variety of image format. For high quality picture, one can increase
% the output resolution by Help -> Preferences -> Resolution. Higher resolution
% will take longer time to copy or save, and result in larger file. If needed,
% one can change to white background for picture output. With white background,
% the threshold for the background image needs to be carefully selected to avoid
% strange effect with some background images. For this reason, white background
% is intended only for picture output.
% 
% The selected NIfTI can also be saved as different format from File -> Save
% NIfTI as. For example, a file can be saved as a different resolution. With a
% transformation matrix, a file can also be saved into a different template. The
% latter is needed for FSLview since it won't allow overlay with different
% resolution or dimension at least till version 5.0.8.
% 
% See also NII_TOOL DICM2NII NII_XFORM

%% By Xiangrui Li (xiangrui.li at gmail.com)
% History(yymmdd):
% 151021 Almost ready to publish.
% 151104 Include drag&drop by Maarten van der Seijs.
% 151105 Bug fix for Show NIfTI hdr/ext.
% 151106 Use hs.q.interp to avoid unnecessary interp for overlay;
%        Use check mark for colorbar/crosshair/background menu items.
% 151111 Take care of slope/inter of img; mask applies to DTI lines.
% 151114 Avoid see-thu for white background.
% 151119 Make smooth faster by using only 3 slices; Allow flip L/R.
% 151120 Implement key navigation and zoom in/out.
% 151121 Implement erode for white background; Add 'Show NIfTI essentials'. 
% 151122 Bug fix for background img with sform=0.
% 151123 Show coordinate system in fig title; show masked/aligned in file list;
%        Bug fix for alpha (0/1 only); Bug fix for background image R change.
% 151125 Avoid recursion for white background (failed under Linux).
% 151128 Change checkbox to radiobutton: looks better;
%        Fix the bug in reorient dim(perm), introduced in last revision.
% 151129 Correct Center crosshair (was off by 0.5); Use evt.Key & add A/C/X/F1. 
% 151201 Keep fig location & pref when 'open' & dnd: more user friendly.
% 151202 java robot trick allows to add overlay by Ctrl-dnd.
% 151206 Bug fix for RGB nii smooth; Change hs.q(i) to hs.q{i}.
% 151207 Implement 'Add modulation'.
% 151208 Add xyz display for each view.
% 151217 Callback uses subfunc directly, and include fh as input.
% 151222 Show ROI labels (eg AAL) if .txt file with same name exists.
% 151224 Implement more matlab LUTs and custom LUT.
% 151229 Store fname in p rather than fh UserData; Remove WheelFcn zoom.
% 151230 Use listbox for files; Add stack buttons; file order reversed.
% 160102 Store hs.q{i}.R0 if need interp; No coordinate change for background.
% 160104 set_cdata: avoid indexing for single vol img: may be much faster!
%        jscroll_handle from findjobj.m to set vertical scoll bar as needed.
% 160107 Rename XYZ label to IJK, implement "Set crosshair at XYZ".
% 160108 Fix the case of 2-form_code for background and addMask.
% 160109 Allow to turn off mouse-over display from Preferences;
%        Implement Help -> Check update for easy update;
%        Use Matlab built-in pref method for all files in the package.
% 160111 Use image ButtonDownFcn; Mouse down&move now is same as mouse move.
% 160112 Change back to line for crosshair: quiver is slow;
%        Bug fix for white backgrnd & saveas backgrnd due to list order reverse.
% 160113 Implement time course for a cube.
% 160114 Figure visible avoids weird hanging problem for some matlab versions.
% 160119 Allow adding multiple overlays with cellstr as 2nd input.
% 160131 set_file: bug fix for cb(j); dnd: restore mouse location.
% 160208 Allow moving background image in stack.
% 160209 RGBA data supported; Background image can use lut "lines".
% 160218 "lines": Fix non-isovoxel display; Same-dim background not required.
% 160402 nii_xform_mat: make up R for possible Analyze file.
% 160506 phase LUT to map complex img: useful for retinotopy.
% 160509 Have 3 phase LUTs; Implement 'Open in new window'.
% 160512 KeyPressFcn bug fix: use the smallest axis dim when zoom in.
% 160517 KeyPressFcn: avoid double-dlg by Ctrl-A; phase6: bug fix for blue b3.
% 160531 use handle() for fh & others: allow dot convention for early matlab.
% 160601 Add cog and to image center, re-organize 'Set crosshair at' menu.
% 160602 bug fix for 'closeAll' files.String(ib); COG uses abs and excludes NaN.
% 160605 Add 'RGB' LUT to force RGB display: DTI_V1 or fsl style RGB file.
% 160608 javaDropFcn: 2 more method for ctlDn; bug fix for fh.xxx in Resize.
% 160620 Use JIDE CheckBoxList; Simplify KeyFcn by not focusing on active items. 
% 160627 javaDropFcn: robot-click drop loc to avoid problem with vnc display. 
% 160701 Implement hist for current volume; hs.value show top overlay first.
% 160703 bug fix for 'Add aligned' complex img: set q.phase after re-orient.
% 160710 Implement 'Create ROI file'; Time coure is for sphere, not cube.
% 160713 javaDropFnc for Linux: Robot key press replaces mouse click;
%        Implement 'Set crosshair at' 'Smoothed maximum'.
% 160715 lut2img: bug fix for custom LUT; custom LUT uses rg [0 255]; add gap=0.
% 160721 Implement 'Set crosshair at' & 'a point with value of'.
% 160730 Allow to load single volume for large dataset.
% 161003 Add aligned overlay: accept FSL warp file as transformation.
% 161010 Implement 'Save volume as'; xyzr2roi: use valid q/sform.
% 161018 Take care of issue converting long file name to var name.
% 161103 Fix qform-only overlay, too long fig title, overlay w/o valid formCode.
% 161108 Implement "Crop below crosshair" to remove excessive neck tissue.
% 161115 Use .mat file for early spm Analyze file.
% 161216 Show more useful 4x4 R for both s/q form in nii essentials.
% 170109 bug fix: add .phase for background nifti.
% 170130 get_range: use nii as input, so always take care of slope/inter.
% 170210 Use flip for flipdim if available (in multiple files).
% 170212 Can open nifti-convertible files; Add Save NIfTI as -> a copy.
% 170421 java_dnd() changed as func, ControlDown OS independent by ACTION_LINK.
% 170515 Use area to normalize histogram.
% 171031 Implement layout. axes replace subplot to avoid overlap problem.
%%

if nargin==2 && ischar(fname) && strcmp(fname, 'func_handle')
    varargout{1} = str2func(overlayName);
    return;
end

pf = getpref('nii_viewer_para');
if isempty(pf) || ~isfield(pf, 'mouseOver')
    pf = struct('openPath', pwd, 'addPath', pwd, 'interp', 'linear', 'extraV', NaN, ...
        'dpi', '0', 'rightOnLeft', false, 'mouseOver', true, 'layout', 2);
    setpref('nii_viewer_para', fieldnames(pf), struct2cell(pf));
end
if ~isfield(pf, 'layout')
    pf.layout = 2;
    setpref('nii_viewer_para', 'layout', 2);
end

if nargin<1
    types = '*.nii; *.hdr; *.nii.gz; *.hdr.gz';
    [fname, pName] = uigetfile([pf.openPath '/' types], ...
        'Select NIfTI to view', 'MultiSelect', 'on');
    if isnumeric(fname), return; end
    fname = strcat([pName '/'], fname);
    setpref('nii_viewer_para', 'openPath', pName);
end

nii = get_nii(fname);
[q, hs.form_code, rg, dim] = read_nii(nii); % re-oriented
q.Ri = inv(q.R);
nVol = size(q.nii.img, 4);
hs.pixdim = q.pixdim;
if ~isreal(q.nii.img)
    q.phase = angle(q.nii.img); % -pi to pi
    q.phase = mod(q.phase/(2*pi), 1); % 0~1
    q.nii.img = abs(q.nii.img); % real now
end
hs.q{1} = q;
hs.dim = single(dim); % single saves memory for ndgrid
[siz, plotPos] = plot_pos(dim.*hs.pixdim, pf.layout);
hs.gap = min(hs.pixdim) ./ hs.pixdim * 3; % in unit of smallest pixdim

p = struct('lut', [], 'lb', rg(1), 'ub', rg(2));
p = dispPara(p, hs.q{1}.nii.hdr);
[pName, niiName, ext] = fileparts(p.fname);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end

res = screen_pixels(1); % use 1st screen
if nargin>1 && all(ishandle(overlayName)) % internal call by 'open' or dnd
    fh = overlayName;
    fn = get(fh, 'Number');
    pos = fh.Position(1:2);
    close(fh); 
else
    pos = round((res-siz)/2);

    set(0, 'ShowHiddenHandles', 'on');
    a = handle(findobj('Type', 'figure', 'Tag', 'nii_viewer'));
    set(0, 'ShowHiddenHandles', 'off');
    if isempty(a)
        fn = 'ni' * 256.^(1:2)'; % start with a big number for figure
    elseif numel(a) == 1
        fn = get(a, 'Number') + 1;
    else
        fn = max(cell2mat(get(a, 'Number'))) + 1;
    end
end
fh = figure(fn);
hs.fig = handle(fh); % have to use numeric for uipanel for older matlab
figNam = p.fname;
if numel(figNam)>40, figNam = [figNam(1:40) '...']; end
figNam = ['nii_viewer - ' figNam ' (' formcode2str(hs.form_code(1)) ')'];
if pos(1)+siz(1) > res(1), pos(1) = res(1)-siz(1)-10; end
if pos(2)+siz(2) > res(2)-130, pos(2) = min(pos(2), 50); end
set(fh, 'Toolbar', 'none', 'Menubar', 'none', 'Renderer', 'opengl', ...
    'NumberTitle', 'off', 'Tag', 'nii_viewer', 'DockControls', 'off', ...
    'Position', [pos siz+[0 64]], 'Name', figNam);
cb = @(cmd) {@nii_viewer_cb cmd hs.fig}; % callback shortcut
xyz = [0 0 0]; % start cursor location
c = round(hs.q{1}.Ri * [xyz 1]'); c = c(1:3)' + 1; % 
ind = c<=1 | c>=dim;
c(ind) = round(dim(ind)/2);
% c = round(dim/2); % start cursor at the center of images
xyz = round(hs.q{1}.R * [c-1 1]'); % take care of rounding error

%% control panel
pos = getpixelposition(fh); pos = [1 pos(4)-64 pos(3) 64];
hs.panel = uipanel(fh, 'Units', 'pixels', 'Position', pos, 'BorderType', 'none');
hs.focus = uicontrol(hs.panel, 'Style', 'text'); % dummy uicontrol for focus

% file list by JIDE CheckBoxList: check/selection independent
mdl = handle(javax.swing.DefaultListModel, 'CallbackProperties'); % dynamic item
mdl.add(0, niiName);
% mdl.IntervalAddedCallback = cb('width'); % seems not needed
mdl.IntervalRemovedCallback = cb('width');
mdl.ContentsChangedCallback = cb('width'); % add '(mask)' etc
h = handle(com.jidesoft.swing.CheckBoxList(mdl), 'CallbackProperties');
h.setFont(java.awt.Font('Tahoma', 0, 11));
% h.ClickInCheckBoxOnly = true; % it's default
h.setSelectionMode(0); % single selection
h.setSelectedIndex(0); % 1st highlighted
h.addCheckBoxListSelectedIndex(0); % check it
h.ValueChangedCallback = {@set_file fh}; % selection change
h.MouseReleasedCallback = @(~,~)uicontrol(hs.focus); % move focus away
h.setToolTipText(['<html>Select image to show/modify its display ' ...
    'parameters.<br>Click checkbox to turn on/off image']);
jScroll = com.mathworks.mwswing.MJScrollPane(h);
width = h.getPreferredScrollableViewportSize.getWidth;
width = max(60, min(width+20, pos(3)-408)); % 20 pixels for vertical scrollbar
[~, hs.scroll] = javacomponent(jScroll, [2 4 width 60], hs.panel);
hCB = handle(h.getCheckBoxListSelectionModel, 'CallbackProperties');
hCB.ValueChangedCallback = cb('toggle'); % check/uncheck
hs.files = javaObjectEDT(h); % trick to avoid java error by Yair

% panel for controls except hs.files
pos = [width 1 pos(3)-width pos(4)];
ph = uipanel(hs.panel, 'Units', 'pixels', 'Position', pos, 'BorderType', 'none');
clr = get(ph, 'BackgroundColor');
hs.params = ph;

feature('DefaultCharacterSet', 'UTF-8'); % for old matlab to show triangles
hs.overlay(1) = uicontrol(ph, 'Style', 'pushbutton', 'FontSize', 7, ...
    'Callback', cb('stack'), 'Enable', 'off', 'SelectionHighlight', 'off', ...
    'String', char(9660), 'Position', [1 36 16 15], 'Tag', 'down', ...
    'TooltipString', 'Move selected image one level down');
hs.overlay(2) = copyobj(hs.overlay(1), ph);
set(hs.overlay(2), 'Callback', cb('stack'), ...
    'String', char(9650), 'Position', [1 50 16 15], 'Tag', 'up', ...
    'TooltipString', 'Move selected image one level up');

hs.value = uicontrol(ph, 'Style', 'text', 'Position', [208 38 pos(3)-208 20], ...
    'BackgroundColor', clr, 'FontSize', 8+(~ispc && ~ismac), ...
    'TooltipString', '(i,j,k)=(x,y,z): top ... bottom');

% IJK java spinners
labls = 'IJK';
str = {'Left to Right' 'Posterior to Anterior' 'Inferior to Superior'};
pos = [38 44 22]; posTxt = [36 10 20];
for i = 1:3
    loc = [(i-1)*64+34 pos];
    txt = sprintf('%s, 1:%g', str{i}, dim(i));
    hs.ijk(i) = java_spinner(loc, [c(i) 1 dim(i) 1], ph, cb('ijk'), '#', txt);
    uicontrol(ph, 'Style', 'text', 'String', labls(i), 'BackgroundColor', clr, ...
        'Position', [loc(1)-11 posTxt], 'TooltipString', txt, 'FontWeight', 'bold');
end

% Controls for each file
uipanel('Parent', ph, 'Units', 'pixels', 'Position', [1 2 408 34], ...
    'BorderType', 'etchedin', 'BorderWidth', 2);
hs.lb = java_spinner([7 8 52 22], [p.lb -inf inf p.lb_step], ph, ...
    cb('lb'), '#.##', 'min value (threshold)');
hs.ub = java_spinner([59 8 52 22], [p.ub -inf inf p.ub_step], ph, ...
    cb('ub'), '#.##', 'max value (clipped)');
hs.lutStr = {'grayscale' 'red' 'green' 'blue' 'violet' 'yellow' 'cyan' ...
    'red-yellow' 'blue-green' 'two-sided'  '<html><font color="red">lines' ...
    'parula' 'jet' 'hsv' 'hot' 'cool' 'spring' 'summer' 'autumn' 'winter' ...
    'bone' 'copper' 'pink' 'prism' 'flag' 'phase' 'phase3' 'phase6' 'RGB'};
hs.lut = uicontrol(ph, 'Style', 'popup', 'Position', [113 8 74 22], ...
    'String', hs.lutStr, 'BackgroundColor', 'w', 'Callback', cb('lut'), ...
    'Value', p.lut, 'TooltipString', 'Lookup table options for non-RGB data');

hs.alpha = java_spinner([187 8 44 22], [1 0 1 0.1], ph, cb('alpha'), '#.#', ...
    'Alpha: 0 transparent, 1 opaque');

hs.smooth = uicontrol(ph, 'Style', 'checkbox', 'Value', p.smooth, ...
    'Position', [231 8 60 22], 'String', 'smooth', ...
    'BackgroundColor', clr, 'Callback', cb('smooth'), ...
    'TooltipString', 'Smooth image in 3D');
hs.interp = uicontrol(ph, 'Style', 'popup', 'Position', [291 8 68 22], ...
    'String', {'nearest' 'linear' 'cubic' 'spline'}, 'Value', p.interp, ...
    'Callback', cb('interp'), 'Enable', 'off', 'BackgroundColor', 'w', ... 
    'TooltipString', 'Interpolation method for overlay');
hs.volume = java_spinner([361 8 44 22], [1 1 nVol 1], ph, cb('volume'), '#', ...
    ['Volume number, 1:' num2str(nVol)]);
hs.volume.setEnabled(nVol>1);

%% Three views: sag, cor, tra
% this panel makes resize easy: axes relative to panel
hs.frame = uipanel(fh, 'Units', 'pixels', 'Position', [1 1 siz], ...
    'BorderType', 'none', 'BackgroundColor', 'k');

for i = 1:3
    j = 1:3; j(j==i) = [];
    hs.ax(i) = axes('Position', plotPos(i,:), 'Parent', hs.frame);
    hs.hsI(i) = image(zeros(dim(j([2 1])), 'single')); hold(hs.ax(i), 'on');
    
    x = [c(j(1))+[-1 1 0 0]*hs.gap(j(1)); 0 dim(j(1))+1 c(j(1))*[1 1]];
    y = [c(j(2))+[0 0 -1 1]*hs.gap(j(2)); c(j(2))*[1 1] 0 dim(j(2))+1];
    hs.cross(i,:) = line(x, y);

    hs.xyz(i) = text(0.02, 0.96, num2str(xyz(i)), 'Units', 'normalized', ...
        'Parent', hs.ax(i), 'FontSize', 12);
end
set(hs.hsI, 'ButtonDownFcn', cb('mousedown'), 'Visible', 'off'); % for copyimg
p.hsI = copyimg(hs);
set(hs.scroll, 'UserData', p); % after copying p.hsI
hs.iback = 1;

labls='ASLSLP'; 
pos = [0.95 0.5; 0.47 0.96; 0 0.5; 0.47 0.96; 0 0.5; 0.47 0.05]; 
for i = 1:numel(labls)
    hs.ras(i) = text(pos(i,1), pos(i,2), labls(i), 'Units', 'normalized', ...
        'Parent', hs.ax(ceil(i/2)), 'FontSize', 12, 'FontWeight', 'bold');
end

% early matlab's colormap works only for axis, so ax(4) is needed.
hs.ax(4) = axes('Position', plotPos(4,:), 'Parent', hs.frame, 'Ylim', [0 1]);
colorbar('Location', 'West', 'Units', 'Normalized');
hs.colorbar = findobj(fh, 'Tag', 'Colorbar'); % trick for early matlab
set(hs.colorbar, 'Visible', 'off', 'UIContextMenu', '', 'EdgeColor', [1 1 1]);
% hs.colorbar = colorbar(hs.ax(4), 'YTicks', [0 0.5 1], 'Color', [1 1 1], ...
%     'Location', 'West', 'PickableParts', 'none', 'Visible', 'off');

% image() reverses YDir. Turn off ax and ticks
set(hs.ax, 'YDir', 'normal', 'Visible', 'off');
set([hs.ras hs.cross(:)' hs.xyz], 'Color', 'b', 'UIContextMenu', '');
try set([hs.ras hs.cross(:)' hs.xyz], 'PickableParts', 'none'); % new matlab
catch, set([hs.ras hs.cross(:)' hs.xyz], 'HitTest', 'off'); % old ones
end

%% menus
h = uimenu(fh, 'Label', '&File');
uimenu(h, 'Label', 'Open', 'Accelerator', 'O', 'UserData', pName, 'Callback', cb('open'));
uimenu(h, 'Label', 'Open in new window', 'Callback', cb('open'));
uimenu(h, 'Label', 'Apply mask', 'Callback', @addMask);
uimenu(h, 'Label', 'Apply modulation', 'Callback', @addMask);
uimenu(h, 'Label', 'Load custom LUT', 'Callback', cb('custom'));
h_savefig = uimenu(h, 'Label', 'Save figure as');
h_saveas = uimenu(h, 'Label', 'Save NIfTI as');
uimenu(h, 'Label', 'Save volume as ...', 'Callback', cb('saveVolume'));
uimenu(h, 'Label', 'Crop below crosshair', 'Callback', cb('cropNeck'));
uimenu(h, 'Label', 'Create ROI file ...', 'Callback', cb('ROI'));
uimenu(h, 'Label', 'Close window', 'Accelerator', 'W', 'Callback', 'close gcf');

uimenu(h_saveas, 'Label', 'SPM 3D NIfTI (one file/pair per volume)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'NIfTI standard RGB (for AFNI, later mricron)', ...
    'Callback', @save_nii_as, 'Separator', 'on');
uimenu(h_saveas, 'Label', 'FSL style RGB (RGB saved in dim 4)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'Old mricron style RGB (RGB saved in dim 3)', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'a copy', 'Callback', @save_nii_as, 'Separator', 'on');
uimenu(h_saveas, 'Label', 'file with a new resolution', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'file matching background', 'Callback', @save_nii_as);
uimenu(h_saveas, 'Label', 'file in aligned template space', 'Callback', @save_nii_as);

fmt = {'pdf' 'eps' 'png' 'jpg' 'tif' 'bmp'};
if ispc, fmt = [fmt 'emf']; end
for i = 1:numel(fmt)
    uimenu(h_savefig, 'Label', fmt{i}, 'Callback', cb('save'));
end

if ispc || ismac
    h = uimenu(fh, 'Label', '&Edit');
    uimenu(h, 'Label', 'Copy figure', 'Separator', 'on', 'Callback', cb('copy'));
end

h_over = uimenu(fh, 'Label', '&Overlay');
hs.add = uimenu(h_over, 'Label', 'Add overlay', 'Accelerator', 'A', ...
    'UserData', pf.addPath, 'Callback', cb('add'));
uimenu(h_over, 'Label', 'Add aligned overlay', 'Callback', cb('add'));

h = uimenu(h_over, 'Label', 'Move selected image', 'Enable', 'off');
uimenu(h, 'Label', 'to top',         'Callback', cb('stack'), 'Tag', 'top');
uimenu(h, 'Label', 'to bottom',      'Callback', cb('stack'), 'Tag', 'bottom');
uimenu(h, 'Label', 'one level up',   'Callback', cb('stack'), 'Tag', 'up');
uimenu(h, 'Label', 'one level down', 'Callback', cb('stack'), 'Tag', 'down');
hs.overlay(3) = h;

hs.overlay(5) = uimenu(h_over, 'Label', 'Remove overlay', 'Accelerator', 'R', ...
    'Callback', cb('close'), 'Enable', 'off');
hs.overlay(4) = uimenu(h_over, 'Label', 'Remove overlays', ...
    'Callback', cb('closeAll'), 'Enable', 'off');

h_view = uimenu(fh, 'Label', '&View');
h = uimenu(h_view, 'Label', 'Zoom in by');
for i = [1 1.2 1.5 2 3 4 5 8 10 20]
    uimenu(h, 'Label', num2str(i), 'Callback', cb('zoom'));
end
h = uimenu(h_view, 'Label', 'Layout', 'UserData', pf.layout);
uimenu(h, 'Label', 'one-row', 'Callback', cb('layout'));
uimenu(h, 'Label', 'two-row sag on right', 'Callback', cb('layout'));
uimenu(h, 'Label', 'two-row sag on left', 'Callback', cb('layout'));
uimenu(h_view, 'Label', 'White background', 'Callback', cb('background'));
hLR = uimenu(h_view, 'Label', 'Right on left side', 'Callback', cb('flipLR'));
uimenu(h_view, 'Label', 'Show colorbar', 'Callback', cb('colorbar'));
uimenu(h_view, 'Label', 'Show crosshair', 'Separator', 'on', ...
    'Checked', 'on', 'Callback', cb('cross'));
h = uimenu(h_view, 'Label', 'Set crosshair at');
uimenu(h, 'Label', 'center of view', 'Callback', cb('viewCenter'));
uimenu(h, 'Label', 'center of image', 'Callback', cb('center'));
uimenu(h, 'Label', 'COG of image', 'Callback', cb('cog'));
uimenu(h, 'Label', 'Smoothed maximum', 'Callback', cb('maximum'));
uimenu(h, 'Label', 'a point [x y z] ...', 'Callback', cb('toXYZ'));
uimenu(h, 'Label', 'a point with value of ...', 'Callback', cb('toValue'));
uimenu(h_view, 'Label', 'Crosshair color', 'Callback', cb('color'));
h = uimenu(h_view, 'Label', 'Crosshair gap');
for i = [0 1 2 3 4 5 6 8 10 20 40]
    str = num2str(i); if i==6, str = [str ' (default)']; end %#ok
    uimenu(h, 'Label', str, 'Callback', cb('gap'));
end
h = uimenu(h_view, 'Label', 'Crosshair thickness');
uimenu(h, 'Label', '0.5 (default)', 'Callback', cb('thickness'));
for i = [0.75 1 2 4 8]
    uimenu(h, 'Label', num2str(i), 'Callback', cb('thickness'));
end

h = uimenu(fh, 'Label', '&Window');
uimenu(h, 'Label', 'Show NIfTI essentials', 'Callback', cb('essential'));
uimenu(h, 'Label', 'Show NIfTI hdr', 'Callback', cb('hdr'));
uimenu(h, 'Label', 'Show NIfTI ext', 'Callback', cb('ext'));
uimenu(h, 'Label', 'DICOM to NIfTI converter', 'Callback', 'dicm2nii', 'Separator', 'on');
th = uimenu(h, 'Label', 'Time course ...', 'Callback', cb('tc'), 'Separator', 'on');
setappdata(th, 'radius', 6);
uimenu(h, 'Label', 'Histogram', 'Callback', cb('hist'));

h = uimenu(fh, 'Label', '&Help');
hs.pref = uimenu(h, 'Label', 'Preferences', 'UserData', pf, ...
    'Callback', @(~,~)pref_dialog(gcbo));
uimenu(h, 'Label', 'Key shortcut', 'Callback', cb('keyHelp'));
uimenu(h, 'Label', 'Show help text', 'Callback', 'doc nii_viewer');
checkUpdate = dicm2nii('', 'checkUpdate', 'func_handle');
uimenu(h, 'Label', 'Check update', 'Callback', @(~,~)checkUpdate('nii_viewer'));
uimenu(h, 'Label', 'About', 'Callback', cb('about'));

%% finalize gui
if isnumeric(fh) % for older matlab
    fh = handle(fh);
    schema.prop(fh, 'Number', 'mxArray'); fh.Number = fn;
    hs.lut = handle(hs.lut);
    hs.frame = handle(hs.frame);
    hs.value = handle(hs.value);
    hs.panel = handle(hs.panel);                   
    hs.params = handle(hs.params);                   
    hs.scroll = handle(hs.scroll);                   
end
guidata(fh, hs); % store handles and data

%% java_dnd based on dndcontrol at matlabcentral/fileexchange/53511
try % panel has JavaFrame in later matlab
    jFrame = handle(hs.frame.JavaFrame.getGUIDEView, 'CallbackProperties');
catch
    warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    jFrame = fh.JavaFrame.getAxisComponent;
end
try java_dnd(jFrame, {@javaDropFcn fh}); catch me, disp(me.message); end

set(fh, 'ResizeFcn', cb('resize'), ... % 'SizeChangedFcn' for later matlab
    'WindowKeyPressFcn', @KeyPressFcn, ...
    'PaperPositionMode', 'auto', 'HandleVisibility', 'Callback');
nii_viewer_cb(fh, [], 'resize', fh); % avoid some weird problem

if pf.mouseOver, set(fh, 'WindowButtonMotionFcn', cb('mousemove')); end
if pf.rightOnLeft, nii_viewer_cb(hLR, [], 'flipLR', fh); end
set_cdata(hs);
set_xyz(hs);

if nargin>1
    if ischar(overlayName)
        addOverlay(overlayName, fh);
    elseif iscellstr(overlayName)
        for i=1:numel(overlayName), addOverlay(overlayName{i}, fh); end
    end
end

if hs.form_code(1)<1
    warndlg(['There is no valid form code in NIfTI. The orientation ' ...
        'labeling is likely meaningless.']);
end

%% callbacks
function nii_viewer_cb(h, ~, cmd, fh)
hs = guidata(fh);
switch cmd
    case 'ijk' % IJK spinner
        ix = find(h == hs.ijk);
        set_cdata(hs, ix);
        set_cross(hs, ix);
        xyz = set_xyz(hs);
        for i = 1:3, set(hs.xyz(i), 'String', xyz(i)); end % need 3 if oblique
    case 'mousedown' % image clicked
        % if ~strcmp(get(fh, 'SelectionType'), 'normal'), return; end
        ax = gca;
        c = get(ax, 'CurrentPoint');
        c = round(c(1, 1:2));
        i = 1:3;
        i(ax==hs.ax(1:3)) = [];
        hs.ijk(i(1)).setValue(c(1));
        hs.ijk(i(2)).setValue(c(2));
    case {'lb' 'ub' 'lut' 'alpha' 'smooth' 'interp' 'volume'}
        if ~strcmp(cmd, 'volume'), uicontrol(hs.focus); end % move away focus
        p = hs.scroll.UserData;
        i = hs.files.getSelectedIndex+1;
        val = get(h, 'Value');
        
        if val == 11 && val~=p(i).lut
            hs.lut.UserData = p(i).lut; % remember old lut
        end
        if strcmp(cmd, 'lut')
            err = false;
            if val == 11 % error check for vector lines
                err = true;
                if size(hs.q{i}.nii.img,4)~=3
                    errordlg('Not valid vector data: dim4 is not 3');
                else
                    a = sum(hs.q{i}.nii.img.^2, 4); a = a(a(:)>1e-4);
                    if any(abs(a-1)>0.1)
                        errordlg('Not valid vector data: squared sum is not 1');
                    else, err = false; % passed all checks
                    end
                end
            elseif any(val == 26:28) % error check for phase img
                err = ~isfield(hs.q{i}, 'phase');
                if err, warndlg('Seleced image is not complex data.'); end
            elseif val == 29 % RGB
                err = size(hs.q{i}.nii.img,4)~=3;
                if err, errordlg('RGB LUT requres 3-volume data.'); end
            end
            if err, hs.lut.Value = p(i).lut; return; end % undo selection
        end
        
        p(i).(cmd) = val;
        hs.scroll.UserData = p;
        if any(strcmp(cmd, {'lut' 'lb' 'ub'})), set_colorbar(hs); end
        if strcmp(cmd, 'volume'), set_xyz(hs); end
        set_cdata(hs);
    case 'resize'
        if isempty(hs), return; end
        cb = fh.ResizeFcn;
        fh.ResizeFcn = []; drawnow; % avoid weird effect
        clnObj = onCleanup(@() set(fh, 'ResizeFcn', cb)); % restore func
        
        posP = getpixelposition(hs.panel); % get old height in pixels
        posF = getpixelposition(fh); % asked position by user
        posI = getpixelposition(hs.frame); % old size
        
        siz = hs.frame.Position(3:4);
        res = screen_pixels(1);
        oldF = round([posI(3) posI(4)+posP(4)]); % previous fig size
        if isequal(oldF, posF(3:4)), return; end % moving without size change
        if all(posF(3:4) >= oldF) % enlarge
            a = max([posF(3) posF(4)-posP(4)] ./ siz) * siz;
            a(1) = min(a(1), res(1)-30); % leave space for MAC dock etc
            a(2) = min(a(2), res(2)-92-posP(4)); % leave space for title bar etc
            a = min(a ./ siz) * siz;
        elseif all(posF(3:4) <= oldF) % shrink
            a = min([posF(3) posF(4)-posP(4)] ./ siz) * siz;
        else % one side enlarge, another side shrink: use old size
            a = posI(3:4);
        end
        d = posF(1)+a(1)-res(1);
        if d>0, posF(1) = posF(1) - d; end
        d = posF(2)+a(2)+posP(4)+92-res(2);
        if d>0, posF(2) = posF(2) - d; end
        posF(1) = max(posF(1), 10);
        posF(2) = max(posF(2), 50);
        posF(3:4) = [a(1) a(2)+posP(4)]; % final figure size
        fh.Position = posF; % done for fig
        
        posP(2) = posF(4)-posP(4)+1; 
        posP(3) = posF(3);
        hs.panel.Position = posP; % done for control panel
        hs.frame.Position = [1 1 a]; % done for image panel
        
        nii_viewer_cb([], [], 'width', fh);
    case 'toggle' % turn on/off NIfTI
        i = h.getAnchorSelectionIndex+1;
        if i<1, return; end
        checked = hs.files.getCheckBoxListSelectedIndices+1;
        p = hs.scroll.UserData(i);
        if p.show == any(checked==i), return; end % no change
        p.show = ~p.show;
        hs.scroll.UserData(i) = p;

        states = {'off' 'on'};
        try %#ok<*TRYNC>
            set(p.hsI, 'Visible', states{p.show+1});
            if p.show, set_cdata(hs); end
            set_xyz(hs);
        end
    case 'mousemove'
        % if ~strcmp(get(fh, 'SelectionType'), 'normal'), return; end
        c = cell2mat(get(hs.ax(1:3), 'CurrentPoint'));
        c = c([1 3 5], 1:2); % 3x2
        x = cell2mat(get(hs.ax(1:3), 'XLim')); 
        y = cell2mat(get(hs.ax(1:3), 'YLim')); 
        I = cell2mat(get(hs.ijk, 'Value'))';
        if     c(1,1)>x(1,1) && c(1,1)<x(1,2) && c(1,2)>y(1,1) && c(1,2)<y(1,2)%sag
            I = [I(1) c(1,:)];
        elseif c(2,1)>x(2,1) && c(2,1)<x(2,2) && c(2,2)>y(2,1) && c(2,2)<y(2,2)%cor
            I = [c(2,1) I(2) c(2,2)];
        elseif c(3,1)>x(3,1) && c(3,1)<x(3,2) && c(3,2)>y(3,1) && c(3,2)<y(3,2)%tra
            I = [c(3,:) I(3)];
        end
        set_xyz(hs, I);
    case 'open' % open on current fig or new fig
        pf = get(hs.pref, 'UserData');
        [fname, pName] = uigetfile([pf.openPath '/*.nii; *.hdr;*.nii.gz; *.hdr.gz'], ...
            'Select a NIfTI to view', 'MultiSelect', 'on');
        if isnumeric(fname), return; end
        fname = strcat([pName '/'], fname);
        if strcmp(get(h, 'Label'), 'Open in new window'), nii_viewer(fname);
        else, nii_viewer(fname, fh);
        end
    case 'add' % add overlay
        pName = get(hs.add, 'UserData');
        label = get(h, 'Label');
        if strcmp(label, 'Add aligned overlay')
            [fname, pName] = uigetfile([pName '/*.nii; *.hdr;*.nii.gz;' ...
                '*.hdr.gz'], 'Select overlay NIfTI');
            if ~ischar(fname), return; end
            fname = fullfile(pName, fname);
            [mtx, pName] = uigetfile([pName '/*.mat;*_warp.nii;*_warp.nii.gz'], ...
                'Select FSL mat file or warp file transforming the nii to background');
            if ~ischar(mtx), return; end
            mtx = fullfile(pName, mtx);
            addOverlay(fname, fh, mtx);
        else
            [fname, pName] = uigetfile([pName '/*.nii; *.hdr;*.nii.gz;' ...
                '*.hdr.gz'], 'Select overlay NIfTI', 'MultiSelect', 'on');
            if ~ischar(fname) && ~iscell(fname), return; end
            nii = get_nii(strcat([pName filesep], fname));
            addOverlay(nii, fh);
        end
        setpref('nii_viewer_para', 'addPath', pName);
    case 'closeAll' % close all overlays
        p = hs.scroll.UserData;
        ind = numel(p):-1:1;
        ind(ind==hs.iback) = [];
        for j = ind
            delete(p(j).hsI); % remove image
            hs.files.getModel.remove(j-1);
        end
        
        hs.scroll.UserData = p(hs.iback);
        hs.iback = 1;
        hs.q(ind) = []; guidata(fh, hs);
        hs.files.setSelectedIndex(0);
        set_xyz(hs);
    case 'close' % close selected overlay
        i = hs.files.getSelectedIndex+1;
        if i==hs.iback, return; end % no touch to background
        p = hs.scroll.UserData;
        delete(p(i).hsI); % 3 view
        if i<hs.iback, hs.iback = hs.iback-1; end
        hs.q(i) = []; guidata(fh, hs);
        p(i) = [];
        hs.scroll.UserData = p;
        
        hs.files.getModel.remove(i-1);
        hs.files.setSelectedIndex(max(0, i-2));
        set_xyz(hs);
    case {'hdr' 'ext' 'essential'} % show hdr ext or essential
        j = hs.files.getSelectedIndex+1;
        if strcmp(cmd, 'hdr')
            hdr = hs.q{j}.nii.hdr;
        elseif strcmp(cmd, 'ext')
            if ~isfield(hs.q{j}.nii, 'ext')
                errordlg('No extension for the selected NIfTI'); 
                return;
            end
            hdr = {};
            for i = 1:numel(hs.q{j}.nii.ext)
                if ~isfield(hs.q{j}.nii.ext(i), 'edata_decoded'), continue; end
                hdr{end+1} = hs.q{j}.nii.ext(i).edata_decoded; %#ok
            end
            if isempty(hdr)
                errordlg('No known extension for the selected NIfTI'); 
                return;
            elseif numel(hdr) == 1, hdr = hdr{1};
            end
        elseif strcmp(cmd, 'essential')
            hdr = nii_essential(hs.q{j});
        end
        nam = hs.files.getModel.get(j-1);
        if ~isstrprop(nam(1), 'alpha'), nam = ['x' nam]; end % like genvarname
        nam(~isstrprop(nam, 'alphanum')) = '_'; % make it valid for var name
        nam = [nam '_' cmd];
        nam = strrep(nam, '__', '_');
        n = numel(nam); nm = namelengthmax;
        if n>nm, nam(nm-4:n-4) = ''; end
        assignin('base', nam, hdr);
        evalin('base', ['openvar ' nam]);
    case 'cross' % show/hide crosshairs and RAS labels
        if strcmp(get(h, 'Checked'), 'on')
            set(h, 'Checked', 'off');
            set([hs.cross(:)' hs.ras hs.xyz], 'Visible', 'off');
        else
            set(h, 'Checked', 'on');
            set([hs.cross(:)' hs.ras hs.xyz], 'Visible', 'on');
        end
    case 'color' % crosshair color
        c = uisetcolor(get(hs.ras(1), 'Color'), 'Pick crosshair color');
        if numel(c) ~= 3, return; end
        set([hs.cross(:)' hs.ras hs.xyz], 'Color', c);
    case 'thickness' % crosshair thickness
        c = strtok(get(h, 'Label'));
        set(hs.cross(:)', 'LineWidth', str2double(c));
    case 'gap' % crosshair gap
        c = str2double(strtok(get(h, 'Label')));
        hs.gap = min(hs.pixdim) ./ hs.pixdim * c / 2;
        guidata(fh, hs);
        set_cross(hs, 1:3);
    case 'copy' % copy figure into clipboard
        set(hs.panel, 'Visible', 'off');
        clnObj = onCleanup(@() set(hs.panel, 'Visible', 'on'));
        pf = get(hs.pref, 'UserData');
        print('-dbitmap', '-noui', ['-r' pf.dpi]);
        % print('-dmeta', '-painters');
    case 'save' % save figure as picture
        ext = get(h, 'Label');
        fmt = ext;
        if strcmp(ext, 'jpg'), fmt = 'jpeg';
        elseif strcmp(ext, 'tif'), fmt = 'tiff';
        elseif strcmp(ext, 'eps'), fmt = 'epsc';
        elseif strcmp(ext, 'emf'), fmt = 'meta';
        end
        [fname, pName] = uiputfile(['*.' ext], 'Input file name to save figure');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        if any(strcmp(ext, {'eps' 'pdf' 'emf'})), render = '-painters';
        else, render = '-opengl';
        end
        pf = get(hs.pref, 'UserData');
        set(hs.panel, 'Visible', 'off');
        clnObj = onCleanup(@() set(hs.panel, 'Visible', 'on'));
        print(fname, render, '-noui', ['-d' fmt], ['-r' pf.dpi], '-cmyk');
    case 'colorbar' % colorbar on/off
        if strcmpi(get(hs.colorbar, 'Visible'), 'on')
            set(hs.colorbar, 'Visible', 'off'); 
            set(h, 'Checked', 'off');
        else
            set(hs.colorbar, 'Visible', 'on'); 
            set(h, 'Checked', 'on');
            set_colorbar(hs);
        end
    case 'about'
        reviseDate = dicm2nii('', 'reviseDate', 'func_handle');
        str = sprintf(['nii_viewer.m by Xiangrui Li\n\n' ...
            'Last updated on 20%s\n\n', ...
            'Feedback to: xiangrui.li@gmail.com\n'], reviseDate(mfilename));
        helpdlg(str, 'About nii_viewer')
    case 'stack'
        uicontrol(hs.focus); % move focus out of buttons
        i = hs.files.getSelectedIndex+1;
        p = hs.scroll.UserData;
        n = numel(p);
        switch get(h, 'Tag') % for both uimenu and pushbutton
            case 'up' % one level up
                if i==1, return; end
                for j = 1:3, uistack(p(i).hsI(j)); end
                ind = [1:i-2 i i-1 i+1:n]; i = i-1;
            case 'down' % one level down
                if i==n, return; end
                for j = 1:3, uistack(p(i).hsI(j), 'down'); end
                ind = [1:i-1 i+1 i i+2:n]; i = i+1;
            case 'top'
                if i==1, return; end
                for j = 1:3, uistack(p(i).hsI(j), 'up', i-1); end
                ind = [i 1:i-1 i+1:n]; i = 1;
            case 'bottom'
                if i==n, return; end
                for j = 1:3, uistack(p(i).hsI(j), 'down', n-i); end
                ind = [1:i-1 i+1:n i]; i = n;
            otherwise
                error('Unknown stack level: %s', get(h, 'Tag'));
        end
        
        str = cell(hs.files.getModel.toArray);
        str = str(ind);
        p = p(ind);
        for j = 1:numel(p), hs.files.getModel.set(j-1, str{j}); end
        hs.files.selectNone; % seems unnecessary
        on = find([p.show]);
        if ~isempty(on), hs.files.setCheckBoxListSelectedIndices(on-1); end
        hs.files.setSelectedIndex(i-1);
        hs.scroll.UserData = p;
        hs.iback = find(hs.iback == ind);
        
        hs.q = hs.q(ind);
        guidata(fh, hs);
        set_xyz(hs);
    case 'zoom'
        m = str2double(get(h, 'Label'));
        a = min(hs.dim) / m;
        if a<1, m = min(hs.dim); end
        set_zoom(m, hs);
    case 'background'
        if strcmp(get(h, 'Checked'), 'on')
            set(h, 'Checked', 'off');
            hs.frame.BackgroundColor = [0 0 0];
            set(hs.colorbar, 'EdgeColor', [1 1 1]);
        else
            set(h, 'Checked', 'on');
            hs.frame.BackgroundColor = [1 1 1];
            set(hs.colorbar, 'EdgeColor', [0 0 0]);
        end
        set_cdata(hs);
    case 'flipLR'
        if strcmp(get(h, 'Checked'), 'on')
            set(h, 'Checked', 'off');
            set(hs.ax([2 3]), 'XDir', 'normal');
            set(hs.ras([3 5]), 'String', 'L');
        else
            set(h, 'Checked', 'on');
            set(hs.ax([2 3]), 'XDir', 'reverse');
            set(hs.ras([3 5]), 'String', 'R');
        end
    case 'layout'
        submenu = {'one-row' 'two-row sag on right' 'two-row sag on left'};
        layout = find(strcmp(get(h, 'Label'), submenu));
        hLayout = findobj(hs.fig, 'Type', 'uimenu', 'Label', 'Layout');
        if get(hLayout, 'UserData') == layout, return; end
        
        cb = hs.fig.ResizeFcn;
        hs.fig.ResizeFcn = ''; drawnow;
        clnObj = onCleanup(@() set(hs.fig, 'ResizeFcn', cb));
        
        res = screen_pixels(1) - [60 92];
        oldSiz = hs.frame.Position(3:4);
        [siz, plotPos] = plot_pos(hs.dim.*hs.pixdim, layout);
        if siz(1)<oldSiz(1), siz = siz/siz(1)*oldSiz(1); end
        if siz(2)+64>res(2), siz = siz/siz(2)*(res(2)-64); end
        pos = hs.fig.Position;
        pos(2) = pos(2) + oldSiz(2) - siz(2);
        pos(3:4) = siz+[0 64];
        d = pos(2)+pos(4) - res(2); 
        if d>0, pos(2) = pos(2) - d; end 
        hs.fig.Position = pos;
        hs.frame.Position(3:4) = siz;
        hs.panel.Position(2) =  hs.fig.Position(4) - 64;
        for i = 1:4, set(hs.ax(i), 'Position', plotPos(i,:)); end
        set(hLayout, 'UserData', layout);
        drawnow;
    case 'keyHelp'
        str = sprintf([ ...
           'Key press available when focus is not in a number dialer:\n\n' ...
           'Left or Right arrow key: Move crosshair left or right.\n\n' ...
           'Up or Down arrow key: Move crosshair superior or inferior.\n\n' ...
           '[ or ] key: Move crosshair posterior or anterior.\n\n' ...
           '< or > key: Decrease or increase volume number.\n\n' ...
           'Ctrl + or - key: Zoom in or out by 10%% around crosshair.\n\n' ...
           'A: Toggle on/off crosshair.\n\n' ...
           'C: Crosshair to view center.\n\n' ...
           'Space: Toggle on/off selected image.\n\n' ...
           'F1: Show help text.\n']);
        helpdlg(str, 'Key Shortcut');
    case 'center' % image center
        jf = hs.files.getSelectedIndex+1;
        dim = hs.q{jf}.nii.hdr.dim(2:4);
        c = round(hs.q{hs.iback}.Ri * (hs.q{jf}.R * [dim/2-1 1]')) + 1;
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'viewCenter'
        c(1) = mean(get(hs.ax(2), 'XLim'));
        c(2) = mean(get(hs.ax(1), 'XLim'));
        c(3) = mean(get(hs.ax(1), 'YLim'));
        c = round(c-0.5);
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'toXYZ'
        c0 = cell2mat(get(hs.ijk, 'Value'));
        c0 = hs.q{hs.iback}.R * [c0-1; 1];
        c0 = sprintf('%g %g %g', round(c0(1:3)));
        str = 'X Y Z coordinates in mm';
        while 1
            a = inputdlg(str, 'Crosshair to xyz', 1, {c0});
            if isempty(a), return; end
            c = sscanf(a{1}, '%g %g %g');
            if numel(c) == 3, break; end
        end
        c = round(hs.q{hs.iback}.Ri * [c(:); 1]) + 1;
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'toValue'
        def = getappdata(h, 'Value');
        if isempty(def), def = 1; end
        def = num2str(def);
        str = 'Input the voxel value';
        while 1
            a = inputdlg(str, 'Crosshair to a value', 1, {def});
            if isempty(a), return; end
            val = sscanf(a{1}, '%g');
            if ~isnan(val), break; end
        end
        setappdata(h, 'Value', val);
        jf = hs.files.getSelectedIndex+1;
        img = hs.q{jf}.nii.img(:,:,:,hs.volume.getValue);
        c = find(img(:)==val, 1);
        if isempty(c)
            nam = strtok(hs.files.getModel.get(jf-1), '(');
            errordlg(sprintf('No value of %g found in %s', val, nam));
            return;
        end
        dim = size(img); dim(numel(dim)+1:3) = 1;
        [c(1), c(2), c(3)] = ind2sub(dim, c); % ijk+1
        c = round(hs.q{hs.iback}.Ri * (hs.q{jf}.R * [c(:)-1; 1])) + 1;
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'cog' % crosshair to img COG
        jf = hs.files.getSelectedIndex+1;
        img = hs.q{jf}.nii.img(:,:,:,hs.volume.getValue);
        c = img_cog(img);
        if any(isnan(c)), errordlg('No valid COG found'); return; end
        c = round(hs.q{hs.iback}.Ri * (hs.q{jf}.R * [c-1; 1])) + 1;
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'maximum' % crosshair to img max
        jf = hs.files.getSelectedIndex+1;
        img = hs.q{jf}.nii.img(:,:,:,hs.volume.getValue);
        img = smooth23(img, 'gaussian', 5);
        img(isnan(img)) = 0;
        img = abs(img);
        [~, I] = max(img(:));
        dim = size(img); dim(end+1:3) = 1;
        c = zeros(3, 1);
        [c(1), c(2), c(3)] = ind2sub(dim, I);
        c = round(hs.q{hs.iback}.Ri * (hs.q{jf}.R * [c-1; 1])) + 1;
        for i = 1:3, hs.ijk(i).setValue(c(i)); end
    case 'custom' % add custom lut
        pName = get(hs.add, 'UserData');
        [fname, pName] = uigetfile([pName '/*.lut'], 'Select a LUT file');
        if ~ischar(fname), return; end
        [~, nam] = fileparts(fname);
        nam = strtok(nam, '.');
        str = hs.lut.String;
        ind = find(strcmp(str, nam));
        if ~isempty(ind)
            set(hs.lut, 'Value', ind);
        else
            fname = fullfile(pName, fname);
            fid = fopen(fname);
            lut = fread(fid, '*uint8');
            fclose(fid);
            lut = reshape(lut, [numel(lut)/3 3]);
            lut = single(lut) / 255;
            hs.lutStr{end+1} = lut;
            guidata(fh, hs);
            set(hs.lut, 'String', [str; nam], 'Value', numel(str)+1);
        end
        nii_viewer_cb(hs.lut, [], 'lut', fh);
    case 'tc' % time course
        i = hs.files.getSelectedIndex+1;
        dim = hs.q{i}.nii.hdr.dim(2:5);
        nam = strtok(hs.files.getModel.get(i-1), '(');
        if dim(4)<2
            errordlg(['There is only 1 volume for ' nam]);
            return;
        end
        
        r = num2str(getappdata(h, 'radius'));
        r = inputdlg('Radius around crosshair (mm):', 'Time course', 1, {r});
        if isempty(r), return; end
        r = str2double(r{1});
        setappdata(h, 'radius', r);
        for j = 3:-1:1, c(j) = get(hs.ijk(j), 'Value'); end % ijk for background
        c = hs.q{hs.iback}.R * [c-1 1]'; % in mm now
        c = c(1:3);
        b = xyzr2roi(c, r, hs.q{i}.nii.hdr); % overlay space
        nv = prod(dim(1:3));
        b = reshape(b, [1 nv]);
        
        img = nii_tool('img', hs.q{i}.nii.hdr.file_name); % may be re-oriented
        dim = size(img);
        img = reshape(img, [nv prod(dim(4:end))])';
        img = img(:, b);
        img = mean(single(img), 2);
        fh1 = figure(mod(fh.Number,10)+i);
        plot(img); 
        xlabel('Volume number');
        c = sprintf('(%g,%g,%g)', round(c));
        set(fh1, 'Name', [nam ' time course around voxel ' c]);
    case 'hist' % plot histgram
        i = hs.files.getSelectedIndex+1;
        if i<1, return; end
        img = hs.q{i}.nii.img(:,:,:,hs.volume.getValue);
        img = sort(img(:));
        img(isnan(img)) = [];
        img(img<hs.lb.getValue) = [];
        img(img>hs.ub.getValue) = [];
        nv = numel(img);
        img0 = unique(img);
        nu = numel(img0);
        n = max([nv/2000 nu/20 10]);
        n = min(round(n), nu);
        if n == nu, edges = img0;
        else, edges = linspace(0,1,n)*double(img(end)-img(1)) + double(img(1));
        end
        nam = strtok(hs.files.getModel.get(i-1), '(');
        fh1 = figure(mod(fh.Number,10)+i);
        set(fh1, 'NumberTitle', 'off', 'Name', nam);
        [y, x] = hist(img, edges);
        bar(x, y/sum(y)/(x(2)-x(1)), 'hist'); % probability density
        xlabel('Voxel values'); ylabel('Probability density');
        title('Histogram between min and max values');
    case 'width' % adjust hs.scroll width
        hs.files.updateUI;
        width = hs.panel.Position(3);
        x = hs.files.getPreferredScrollableViewportSize.getWidth;
        x = max(60, min(x+20, width-408)); % 408 width of the little panel
        hs.scroll.Position(3) = x;
        hs.params.Position(1) = x+2;
        hs.value.Position(3) = max(1, width-x-hs.value.Position(1));
    case 'saveVolume' % save 1 or more volumes as a nifti
        jf = hs.files.getSelectedIndex+1;
        nam = hs.scroll.UserData(jf).fname;
        t = hs.scroll.UserData(jf).volume;
        while 1
            a = inputdlg('Volume indice to save (2:4 for example)', ...
                'Save Volume', 1, {num2str(t)});
            if isempty(a), return; end
            try
                t = eval(['[' a{1} '];']);
                break;
            end
        end
        pName = fileparts(nam);
        [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
            'Input name to save volume as');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        
        nii = nii_tool('load', nam); % re-load to be safe
        nii.img =  nii.img(:,:,:,t);
        nii_tool('save', nii, fname);
    case 'ROI' % save sphere
        c0 = cell2mat(get(hs.ijk, 'Value'));
        c0 = hs.q{hs.iback}.R * [c0-1; 1];
        c0 = sprintf('%g %g %g', round(c0(1:3)));
        str = {'X Y Z coordinates in mm' 'Radius in mm'};
        while 1
            a = inputdlg(str, 'Sphere ROI', 1, {c0 '6'});
            if isempty(a), return; end
            c = sscanf(a{1}, '%g %g %g');
            r = sscanf(a{2}, '%g');
            if numel(c) == 3, break; end
        end
        
        pName = get(hs.add, 'UserData');
        [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
            'Input file name to save ROI into');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        
        nii = hs.q{hs.files.getSelectedIndex+1}.nii; % original hdr
        b = xyzr2roi(c, r, nii.hdr);        
        nii.img = single(b); % single better supported by FSL
        nii_tool('save', nii, fname);
    case 'cropNeck'
        k0 = get(hs.ijk(3), 'Value') - 1;
        jf = hs.files.getSelectedIndex+1;
        nam = hs.scroll.UserData(jf).fname;
        pName = fileparts(nam);
        [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
            'Input file name to save cropped image');
        if ~ischar(fname), return; end
        fname = fullfile(pName, fname);
        
        d = single(hs.q{jf}.nii.hdr.dim(2:4));
        I = ones([d 4], 'single');
        [I(:,:,:,1), I(:,:,:,2), I(:,:,:,3)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
        I = permute(I, [4 1 2 3]);
        I = reshape(I, [4 prod(d)]); % ijk in 4 by nVox
        R = nii_xform_mat(hs.q{jf}.nii.hdr, hs.form_code); % original R
        k = hs.q{hs.iback}.Ri * R * I; % background ijk
        
        nii = nii_tool('load', nam);
        d = size(nii.img);
        n = prod(d(1:3));
        nii.img = reshape(nii.img, [n prod(d)/n]);
        nii.img(k(3,:)<k0, :) = 0;
        nii.img = reshape(nii.img, d);
        nii_tool('save', nii, fname);
    otherwise
        error('Unknown Callback: %s', cmd);
end

%% zoom in/out with a factor
function set_zoom(m, hs)
c = hs.dim(:) / 2;
if m <= 1, I = c; % full view regardless of crosshair location
else, I = cell2mat(get(hs.ijk, 'Value'));
end
lim = round([I I] + c/m*[-1 1]) + 0.5;
axis(hs.ax(1), [lim(2,:) lim(3,:)]);
axis(hs.ax(2), [lim(1,:) lim(3,:)]);
axis(hs.ax(3), [lim(1,:) lim(2,:)]);

%% KeyPressFcn for figure
function KeyPressFcn(fh, evt)
set(fh, 'UserData', evt.Modifier); % for javaDropFnc
if any(strcmp(evt.Key, evt.Modifier)), return; end % only modifier
hs = guidata(fh);
if ~isempty(intersect({'control' 'command'}, evt.Modifier))
    switch evt.Key
        case {'add' 'equal'}
            [dim, i] = min(hs.dim);
            if     i==1, d = get(hs.ax(2), 'XLim');
            elseif i==2, d = get(hs.ax(1), 'XLim');
            else,        d = get(hs.ax(1), 'YLim');
            end
            d = abs(diff(d'));
            if d<=3, return; end % ignore
            m = dim / d * 1.1;
            if round(dim/2/m)==d/2, m = dim / (d-1); end
            set_zoom(m, hs);
        case {'subtract' 'hyphen'}
            d = abs(diff(get(hs.ax(2), 'XLim')));
            m = hs.dim(1) / d;
            if m<=1, return; end
            m = m / 1.1;
            if round(hs.dim(1)/2/m)==d/2, m = hs.dim(1) / (d+1); end
            if m<1.01, m = 1; end
            set_zoom(m, hs);
    end
    return;
end

switch evt.Key
    case 'leftarrow'
        val = max(get(hs.ijk(1), 'Value')-1, 1);
        hs.ijk(1).setValue(val);
    case 'rightarrow'
        val = min(get(hs.ijk(1), 'Value')+1, hs.dim(1));
        hs.ijk(1).setValue(val);
    case 'uparrow'
        val = min(get(hs.ijk(3),'Value')+1, hs.dim(3));
        hs.ijk(3).setValue(val);
    case 'downarrow'
        val = max(get(hs.ijk(3),'Value')-1, 1);
        hs.ijk(3).setValue(val);
    case 'rightbracket' % ]
        val = min(get(hs.ijk(2),'Value')+1, hs.dim(2));
        hs.ijk(2).setValue(val);
    case 'leftbracket' % [
        val = max(get(hs.ijk(2),'Value')-1, 1);
        hs.ijk(2).setValue(val);
    case 'period' % . or >
        val = min(get(hs.volume,'Value')+1, get(hs.volume.Model,'Maximum'));
        hs.volume.setValue(val);
    case 'comma' % , or <
        val = max(get(hs.volume,'Value')-1, 1);
        hs.volume.setValue(val);
    case 'c'
        nii_viewer_cb([], [], 'viewCenter', hs.fig);
    case {'x' 'space'}
        i = hs.files.getSelectedIndex;
        checked = hs.files.getCheckBoxListSelectedIndices;
        if any(i == checked)
            hs.files.removeCheckBoxListSelectedIndex(i);
        else
            hs.files.addCheckBoxListSelectedIndex(i);
        end
    case 'a'
        h = findobj(hs.fig, 'Type', 'uimenu', 'Label', 'Show crosshair');
        nii_viewer_cb(h, [], 'cross', hs.fig);
    case 'f1'
        doc nii_viewer;
    case 'tab' % prevent tab from cycling uicontrol
        mousexy = get(0, 'PointerLocation'); % for later restore
        posF = getpixelposition(fh);
        posA = getpixelposition(hs.ax(4), true); % relative to figure
        c = posF(1:2) + posA(1:2) + posA(3:4)/2; % ax(4) center xy
        res = screen_pixels;
        rob = java.awt.Robot();
        rob.mouseMove(c(1), res(2)-c(2));
        rob.mousePress(16); rob.mouseRelease(16); % BUTTON1
        set(0, 'PointerLocation', mousexy); % restore mouse location
end

%% Drop callback: drop as background, Ctrl-drop as overlay
function javaDropFcn(~, evt, fh)
try
    nii = get_nii(evt.Data);
    if evt.ControlDown, addOverlay(nii, fh);
    else, nii_viewer(nii, fh);
    end
catch me
    errordlg(me.message);
end

%% update CData/AlphaData for 1 or 3 of the sag/cor/tra views
function set_cdata(hs, iaxis)
if nargin<2, iaxis = 1:3; end
interStr = get(hs.interp, 'String');
p = hs.scroll.UserData;
for i = 1:numel(p)
    if ~p(i).show, continue; end % save time, but need to update when enabled
    lut = p(i).lut;
    if lut == 11 % "lines" special case: do it separately
        vector_lines(hs, i, iaxis); continue; 
    elseif ~strcmpi(get(p(i).hsI(1), 'Type'), 'image') % was "lines"
        delete(p(i).hsI); % delete quiver
        p(i).hsI = copyimg(hs);
        hs.scroll.UserData = p; % update hsI
        if i>1, for j=1:3; uistack(p(i).hsI(j), 'down', i-1); end; end
    end
    t = round(p(i).volume);
    img = hs.q{i}.nii.img;
    isRGB = size(img, 8)>2;
    if isRGB % avoid indexing for single vol img: could speed up a lot
        img = permute(img(:,:,:,t,:,:,:,:), [1:3 8 4:7]);
    elseif size(img,4)>1 && lut~=29
        img = img(:,:,:,t);
    end    
    if ~isfloat(img)
        img = single(img);
        if isfield(hs.q{i}, 'scl_slope')
            img = img * hs.q{i}.scl_slope + hs.q{i}.scl_inter;
        end
    end
    
    if isfield(hs.q{i}, 'mask')
        img = bsxfun(@times, img, hs.q{i}.mask);
    end
    if isfield(hs.q{i}, 'modulation')
        img = bsxfun(@times, img, hs.q{i}.modulation);
    end

    if any(lut == 26:28) % interp/smooth both mag and phase
        img(:,:,:,2) = hs.q{i}.phase(:,:,:,t);
    end
    
    dim4 = size(img,4);
    for ix = iaxis
        ind = get(hs.ijk(ix), 'Value'); % faster than hs.ijk(ix).getValue
        ind = round(ind);
        if ind<1 || ind>hs.dim(ix), continue; end
        ii = {':' ':' ':'};
        io = ii;
        d = hs.dim;
        d(ix) = 1; % 1 slice at dim ix
        im = zeros([d dim4], 'single');
        
        if isfield(hs.q{i}, 'R0') % interp, maybe smooth
            I = ones([d 4], 'single');
            [I(:,:,:,1), I(:,:,:,2), I(:,:,:,3)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
            I = permute(I, [4 1 2 3]);
            I = reshape(I, [4 prod(d)]); % ijk grids of background img
            I(ix,:) = ind-1;
            
            if isfield(hs.q{i}, 'warp')
                iw = {':' ':' ':' ':'}; iw{ix} = ind;
                warp = hs.q{i}.warp(iw{:});
                warp = reshape(warp, [prod(d) 3])'; warp(4,:) = 0;
                I = hs.q{i}.Ri * (hs.q{i}.R0 * I + warp) + 1;
            else
                I = hs.q{i}.Ri * hs.q{i}.R0 * I + 1; % ijk+1 for overlay img
            end
            
            for j = 1:dim4
                if p(i).smooth
                    ns = 3; % number of voxels (odd) used for smooth
                    d3 = d; d3(ix) = ns; % e.g. 3 slices
                    b = zeros(d3, 'single');
                    I0 = I(1:3,:);
                    for k = 1:ns % interp for each slice
                        I0(ix,:) = I(ix,:) - (ns+1)/2 + k;
                        a = interp3a(img(:,:,:,j), I0, interStr{p(i).interp});
                        ii{ix} = k; b(ii{:}) = reshape(a, d);
                    end
                    b = smooth23(b, 'gaussian', ns);
                    io{ix} = (ns+1)/2;
                    im(:,:,:,j) = b(io{:}); % middle one
                else
                    a = interp3a(img(:,:,:,j), I, interStr{p(i).interp});
                    im(:,:,:,j) = reshape(a, d);
                end
            end
        elseif p(i).smooth % smooth only
            ns = 3; % odd number of slices to smooth
            ind1 = ind - (ns+1)/2 + (1:ns); % like ind+(-1:1)
            if any(ind1<1 | ind1>hs.dim(ix)), ind1 = ind; end % 2D
            ii{ix} = ind1;
            io{ix} = mean(1:numel(ind1)); % middle slice
            for j = 1:dim4
                a = smooth23(img(ii{:},j), 'gaussian', ns);
                im(:,:,:,j) = a(io{:});
            end
        else % no interp or smooth
            io{ix} = ind;
            im(:) = img(io{:}, :);
        end
        
        if     ix == 1, im = permute(im, [3 2 4 1]);
        elseif ix == 2, im = permute(im, [3 1 4 2]);
        elseif ix == 3, im = permute(im, [2 1 4 3]);
        end
        
        if ~isRGB % not NIfTI RGB
            rg = sort([p(i).lb p(i).ub]);
            if lut==29 && hs.q{i}.nii.hdr.datatype==2, rg = [0 255]; end % uint8
            [im, alfa] = lut2img(im, lut, rg, hs.lutStr{lut});
        elseif dim4 == 3 % NIfTI RGB
            if max(im(:))>2, im = im / 255; end % guess uint8
            im(im>1) = 1; im(im<0) = 0;
            alfa = sum(im,3) / dim4; % avoid mean
        elseif dim4 == 4 % NIfTI RGBA
            if max(im(:))>2, im = im / 255; end % guess uint8
            im(im>1) = 1; im(im<0) = 0;
            alfa = im(:,:,4);
            im = im(:,:,1:3);
        else
            error('Unknown data type: %s', p(i).fname);
        end
        
        if i==hs.iback && isequal(hs.frame.BackgroundColor, [1 1 1])
            alfa = img2mask(alfa);
        elseif dim4 ~= 4
            alfa = alfa > 0;
        end
        alfa = p(i).alpha * alfa;
        set(p(i).hsI(ix), 'CData', im, 'AlphaData', alfa);
    end
end

%% Add an overlay
function addOverlay(fname, fh, mtx)
hs = guidata(fh);
frm = hs.form_code;
aligned = nargin>2;
R_back = hs.q{hs.iback}.R;
if ~exist('flip', 'builtin'), eval('flip=@flipdim;'); end
if aligned % aligned mtx: do it in special way
    [q, ~, rg, dim] = read_nii(fname, frm, 0); % no re-orient
    R0 = nii_xform_mat(hs.q{hs.iback}.nii.hdr, frm(1)); % original background R
    
    try
        [~, ~, ext] = fileparts(mtx);
        if strcmpi(ext, '.mat')
            R = load(mtx, '-ascii');
            if ~isequal(size(R), [4 4])
                error('Invalid transformation matrix file: %s', mtx);
            end
        else % see nii_xform
            R = eye(4);
            warp = nii_tool('img', mtx);
            if ~isequal(size(warp), [hs.q{hs.iback}.nii.hdr.dim(2:4) 3])
                error('warp file and template file img size don''t match.');
            end
            if det(R0(1:3,1:3))<0, warp(:,:,:,1) = -warp(:,:,:,1); end
            [~, perm, flp] = reorient(R0, hs.q{hs.iback}.nii.hdr.dim(2:4));
            if ~isequal(perm, 1:3), warp = permute(warp, [perm 4]); end
            for j = 1:3
                if flp(j), warp = flip(warp, j); end
            end
            q.warp = warp;
            q.R0 = R_back; % always interp
        end
    catch me
        errordlg(me.message);
        return;
    end

    % see nii_xform for more comment on following method
    R = R0 / diag([hs.q{hs.iback}.nii.hdr.pixdim(2:4) 1]) * R * diag([q.pixdim 1]);
    [~, i1] = max(abs(q.R(1:3,1:3)));
    [~, i0] = max(abs(R(1:3,1:3)));
    flp = sign(R(i0+[0 4 8])) ~= sign(q.R(i1+[0 4 8]));
    if any(flp)
        rotM = diag([1-flp*2 1]);
        rotM(1:3,4) = (dim-1).* flp;
        R = R / rotM;
    end
            
    [q.R, perm, q.flip] = reorient(R, dim); % in case we apply mask to it
    if ~isequal(perm, 1:3)
        dim = dim(perm);
        q.pixdim = q.pixdim(perm);
        q.nii.img = permute(q.nii.img, [perm 4:8]);
    end
    for j = 1:3
        if q.flip(j), q.nii.img = flip(q.nii.img, j); end
    end
    q.alignMtx = mtx; % info only for NIfTI essentials
else
    [q, frm, rg, dim] = read_nii(fname, frm);

    % Silently use another background R_back matching overlay: very rare
    if frm>0 && frm ~= hs.form_code(1) && any(frm == hs.form_code)
        R = nii_xform_mat(hs.q{hs.iback}.nii.hdr, frm); % img alreay re-oriented
        R_back = reorient(R, hs.q{hs.iback}.nii.hdr.dim(2:4));
    elseif frm==0 && isequal(q.nii.hdr.dim(2:4), hs.q{hs.iback}.nii.hdr.dim(2:4))
        q.R = hs.q{hs.iback}.R; % best guess: use background xform
        q.perm = hs.q{hs.iback}.perm;
        q.nii.img = permute(q.nii.img, [q.perm 4:8]);
        for i = 1:3
            if q.flip(i) ~= hs.q{hs.iback}.flip(i)
                q.nii.img = flip(q.nii.img, i);
            end
        end
        q.flip = hs.q{hs.iback}.flip;
        warndlg(['There is no valid coordinate system for the overlay. ' ...
         'The viewer supposes the same coordinate as the background.'], ...
         'Transform Inconsistent');
    elseif frm ~= 2 && ~any(frm == hs.form_code)
        warndlg(['There is no matching coordinate system between the overlay ' ...
         'image and the background image. The overlay is likely meaningless.'], ...
         'Transform Inconsistent');
    end
end

singleVol = 0;
nv = size(q.nii.img, 4);
if nv>1 && numel(q.nii.img)>1e7 % load all or a single volume
    str = ['Input ''all'' or a number from 1 to ' num2str(nv)];
    while 1
        a = inputdlg(str, 'Volumes to load', 1, {'all'});
        if isempty(a), return; end
        a = strtrim(a{1});
        if ~isstrprop(a, 'digit'), break; end
        a = str2num(a);
        if isequal(a,1:nv) || (numel(a)==1 && a>=1 && a<=nv && mod(a,1)==0)
            break; 
        end
    end
    if isnumeric(a) && numel(a)==1
    	singleVol = a;
        q.nii.img = q.nii.img(:,:,:,a);
        rg = get_range(q.nii);
    end
end

ii = [1 6 11 13:15]; % diag and offset: avoid large ratio due to small value
if ~isequal(hs.dim, dim) || any(abs(R_back(ii)./q.R(ii)-1) > 0.01)
    q.R0 = R_back;
end
q.Ri = inv(q.R);

if ~isreal(q.nii.img)
    q.phase = angle(q.nii.img); % -pi to pi
    q.phase = mod(q.phase/(2*pi), 1); % 0~1
    q.nii.img = abs(q.nii.img); % real now
end

% add a slot, and make first for new overlay after no error to read file
p = hs.scroll.UserData;
n = numel(p);
p(2:n+1) = p(1:n);
hs.q = [{q} hs.q];

p(1).hsI = copyimg(hs); % duplicate image obj for overlay: will be at top
p(1).lb = rg(1); p(1).ub = rg(2);
p = dispPara(p, q.nii.hdr);

[pName, niiName, ext] = fileparts(p(1).fname);
if strcmpi(ext, '.gz'), [~, niiName] = fileparts(niiName); end
if aligned, niiName = [niiName '(aligned)']; end
if singleVol, niiName = [niiName '(' num2str(singleVol) ')']; end

try
    mdl = hs.files.getModel;
    mdl.add(n, mdl.get(n-1)); % move to last
    for i = n-1:-1:1, mdl.set(i, mdl.get(i-1)); end % shift by 1
    mdl.set(0, niiName); % new overlay at top
    checked = find([p.show]) - 1;
    hs.files.setCheckBoxListSelectedIndices(checked);
    hs.files.setSelectedIndex(0); hs.files.updateUI;
catch me
    delete(p(1).hsI); % no mess-up in case of error
    errordlg(me.message);
    return;
end
hs.scroll.UserData = p;
set(hs.add, 'UserData', pName);
hs.iback = hs.iback + 1;
guidata(fh, hs); % update hs.q for nii

set_file([], [], fh);
set_cdata(hs);
set_xyz(hs);

%% Reorient 4x4 R
function [R, perm, flp] = reorient(R, dim, leftHand)
% [R, perm, flip] = reorient(R, dim, leftHand)
% Re-orient transformation matrix R (4x4), so it will be diagonal major and
% positive at diagonal, unless the optional third input is true, which requires
% left-handed matrix, where R(1,1) will be negative. 
% The second input is the img space dimension (1x3). 
% The perm output, like [1 2 3] or a permutation of it, indicates if input R was
% permuted for 3 axis. The third output, flip (1x3 logical), indicates an axis 
% (AFTER perm) is flipped if true.
a = abs(R(1:3,1:3));
[~, ixyz] = max(a);
if ixyz(2) == ixyz(1), a(ixyz(2),2) = 0; [~, ixyz(2)] = max(a(:,2)); end
if any(ixyz(3) == ixyz(1:2)), ixyz(3) = setdiff(1:3, ixyz(1:2)); end
[~, perm] = sort(ixyz);
R(:,1:3) = R(:,perm);
flp = R([1 6 11]) < 0; % diag(R(1:3, 1:3))
if nargin>2 && leftHand, flp(1) = ~flp(1); end
rotM = diag([1-flp*2 1]);
rotM(1:3, 4) = (dim(perm)-1) .* flp; % 0 or dim-1
R = R / rotM; % xform matrix after flip

%% Load, re-orient nii, extract essential nii stuff
% nii.img may be re-oriented, but nii.hdr is not touched
function [q, frm, rg, dim] = read_nii(fname, ask_code, reOri)
if ischar(fname), q.nii = nii_tool('load', fname);
else, q.nii = fname; fname = q.nii.hdr.file_name;
end
ndim = q.nii.hdr.dim(1);
dim = q.nii.hdr.dim(2:8);
dim(dim<1 | dim>32767 | mod(dim,1)>0) = 1;
if ndim>4 % 4+ dim, put all into dim4
    if sum(dim(4:7)>1)>1
        warndlg([fname ' has 5 or more dimension. Dimension above 4 are ' ...
            'all treated as volumes for visualization']);        
    end
    dim(4) = prod(dim(4:7)); dim(5:7) = 1;
    q.nii.img = reshape(q.nii.img, [dim size(q.nii.img, 8)]);
end

if nargin<2, ask_code = []; end
[q.R, frm] = nii_xform_mat(q.nii.hdr, ask_code);
dim = dim(1:3);
q.pixdim = q.nii.hdr.pixdim(2:4);
if nargin<3 || reOri
    [q.R, q.perm, q.flip] = reorient(q.R, dim);
    if ~isequal(q.perm, 1:3)
        dim = dim(q.perm);
        q.pixdim = q.pixdim(q.perm);
        q.nii.img = permute(q.nii.img, [q.perm 4:8]);
    end
    if ~exist('flip', 'builtin'), eval('flip=@flipdim;'); end
    for i = 1:3, if q.flip(i), q.nii.img = flip(q.nii.img, i); end; end
else
    q.perm = 1:3;
    q.flip = false(1,3);
end

if size(q.nii.img,4)<4 && ~isfloat(q.nii.img)
    q.nii.img = single(q.nii.img); 
end
if q.nii.hdr.scl_slope==0, q.nii.hdr.scl_slope = 1; end
if q.nii.hdr.scl_slope~=1 || q.nii.hdr.scl_inter~=0
    if isfloat(q.nii.img)
        q.nii.img = q.nii.img * q.nii.hdr.scl_slope + q.nii.hdr.scl_inter;
    else
        q.scl_slope = q.nii.hdr.scl_slope;
        q.scl_inter = q.nii.hdr.scl_inter;
    end
end

% check if ROI labels available: the same file name with .txt extension
if q.nii.hdr.intent_code == 1002 % Label
    [pth, nam, ext] = fileparts(q.nii.hdr.file_name);
    nam1 = fullfile(pth, [nam '.txt']);
    if strcmpi(ext, '.gz') && ~exist(nam1, 'file')
        [~, nam] = fileparts(nam);
        nam1 = fullfile(pth, [nam '.txt']);
    end
    if exist(nam1, 'file') % each line format: 1 ROI_1
        fid = fopen(nam1);
        while 1
            ln = fgetl(fid);
            if ~ischar(ln), break; end
            [ind, a] = strtok(ln);
            ind = str2double(ind);
            try q.labels{ind} = strtrim(a); catch, end
        end
        fclose(fid);
    end
end
rg = get_range(q.nii, isfield(q, 'labels'));

%% Return xform mat and form_code: form_code may have two if not to ask_code
function [R, frm] = nii_xform_mat(hdr, ask_code)
% [R, form] = nii_xform_mat(hdr, asked_code);
% Return the transformation matrix from a NIfTI hdr. By default, this returns
% the sform if available. If the optional second input, required form code, is
% provided, this will try to return matrix for that form code. The second
% optional output is the form code of the actually returned matrix.
fs = [hdr.sform_code hdr.qform_code]; % sform preferred
if fs(1)==fs(2), fs = fs(1); end % sform if both are the same
f = fs(fs>=1 & fs<=4); % 1/2/3/4 only
if isempty(f) || ~strncmp(hdr.magic, 'n', 1) % treat it as Analyze
    frm = 0;
    try % try spm style Analyze
        [pth, nam, ext] = fileparts(hdr.file_name);
        if strcmpi(ext, '.gz'), [~, nam] = fileparts(nam); end
        R = load(fullfile(pth, [nam '.mat']));
        R = R.M;
    catch % make up R for Analyze: suppose xyz order with left storage 
        R = [diag(hdr.pixdim(2:4)) -(hdr.dim(2:4).*hdr.pixdim(2:4)/2)'; 0 0 0 1];
        R(1,:) = -R(1,:); % use left handed
    end
    return;
end

if numel(f)==1 || nargin<2 || isempty(ask_code) % only 1 avail or no ask_code
    frm = f;
else % numel(f) is 2, numel(ask_code) can be 1 or 2
    frm = f(f == ask_code(1));
    if isempty(frm) && numel(ask_code)>1, frm = f(f == ask_code(2)); end
    if isempty(frm) && any(f==2), frm = 2; end % use confusing code 2
    if isempty(frm), frm = f(1); end % no match to ask_code, use sform
end

if frm(1) == fs(1) % match sform_code or no match
    R = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
else % match qform_code
    R = quat2R(hdr);
end

%%
function R = quat2R(hdr)
% Return 4x4 qform transformation matrix from nii hdr.
b = hdr.quatern_b;
c = hdr.quatern_c;
d = hdr.quatern_d;
a = sqrt(1-b*b-c*c-d*d);
if ~isreal(a), a = 0; end % avoid complex due to precision
R = [1-2*(c*c+d*d)  2*(b*c-d*a)     2*(b*d+c*a);
     2*(b*c+d*a)    1-2*(b*b+d*d)   2*(c*d-b*a);
     2*(b*d-c*a )   2*(c*d+b*a)     1-2*(b*b+c*c)];
if hdr.pixdim(1)<0, R(:,3) = -R(:,3); end % qfac
R = R * diag(hdr.pixdim(2:4));
R = [R [hdr.qoffset_x hdr.qoffset_y hdr.qoffset_z]'; 0 0 0 1];

%% Create java SpinnerNumber
function h = java_spinner(pos, val, parent, callback, fmt, helpTxt)
% h = java_spinner(pos, val, parent, callback, fmt, helpTxt)
%  pos: [left bottom width height]
%  val: [curVal min max step]
%  parent: figure or panel
%  fmt: '#' for integer, or '#.#', '#.##'
mdl = javax.swing.SpinnerNumberModel(val(1), val(2), val(3), val(4));
% jSpinner = javax.swing.JSpinner(mdl);
jSpinner = com.mathworks.mwswing.MJSpinner(mdl);
h = javacomponent(jSpinner, pos, parent);
set(h, 'StateChangedCallback', callback, 'ToolTipText', helpTxt);
jEditor = javaObject('javax.swing.JSpinner$NumberEditor', h, fmt);
h.setEditor(jEditor);
h.setFont(java.awt.Font('Tahoma', 0, 11));

%% Estimate lower and upper bound of img display
function rg = get_range(nii, isLabel)
if size(nii.img, 8)>2 || any(nii.hdr.datatype == [128 511 2304]) % RGB / RGBA
    if max(nii.img(:))>2, rg = [0 255]; else, rg = [0 1]; end
    return;
end

img = nii.img(:,:,:,1);
img = img(:);
img(isnan(img) | isinf(img)) = [];
if ~isreal(img), img = abs(img); end
if ~isfloat(img)
    slope = nii.hdr.scl_slope; if slope==0, slope = 1; end
    img = single(img) * slope + nii.hdr.scl_inter;
end

mi = min(img); ma = max(img);
if nii.hdr.intent_code > 1000 || (nargin>1 && isLabel)
    rg = [mi ma]; return;
end

ind = abs(img)>50;
if sum(ind)<numel(img)/10, ind = abs(img)>std(img)/2; end
im = img(ind);
mu = mean(im);
sd = std(im);
rg = mu + [-2 2]*sd;
if rg(1)<=0 && mu-sd>0, rg(1) = sd/5; end
if rg(1)<mi || isnan(rg(1)), rg(1) = mi; end
if rg(2)>ma || isnan(rg(2)), rg(2) = ma; end
if rg(1)==rg(2), rg(1) = mi; end
% rg = round(rg, 2, 'significant'); % since 2014b
rg = str2num(sprintf('%.2g ', rg)); %#ok<*ST2NM>
if rg(1)==rg(2), rg(1) = mi; end
if abs(rg(1))>10, rg(1) = floor(rg(1)/2)*2; end % even number
if abs(rg(2))>10, rg(2) = ceil(rg(2)/2)*2; end % even number

%% Draw vector lines, called by set_cdata
function vector_lines(hs, i, iaxis)
p = hs.scroll.UserData;
d = single(size(hs.q{i}.nii.img));
if strcmpi(get(p(i).hsI(1), 'Type'), 'image') % just switched to "lines"
    delete(p(i).hsI);
    lut = hs.lut.UserData; % last lut
    if isempty(lut), lut = 2; end % default red
    clr = lut2map(lut, hs); clr = clr(end,:);
    cb = get(hs.hsI(1), 'ButtonDownFcn');
    for j = 1:3
        p(i).hsI(j) = quiver(hs.ax(j), 1, 1, 0, 0, 'Color', clr, ...
            'ShowArrowHead', 'off', 'AutoScale', 'off', 'ButtonDownFcn', cb);
    end
    crossFront(hs); % to be safe before next
    if i>1, for j = 1:3, uistack(p(i).hsI(j), 'down', i-1); end; end
    hs.scroll.UserData = p;
    
    if isfield(hs.q{i}, 'R0') && ~isfield(hs.q{i}, 'ivec')
        I = ones([d(1:3) 4], 'single');
        [I(:,:,:,1), I(:,:,:,2), I(:,:,:,3)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
        I = permute(I, [4 1 2 3]);
        I = reshape(I, [4 prod(d(1:3))]);
        I = hs.q{i}.R0 \ (hs.q{i}.R * I) + 1;
        hs.q{i}.ivec = reshape(I(1:3,:)', d);

        R0 = hs.q{i}.R0(1:3, 1:3);
        R0 = bsxfun(@rdivide, R0, sqrt(sum(R0.^2)));
        R = hs.q{i}.R(1:3, 1:3);
        R = bsxfun(@rdivide, R, sqrt(sum(R.^2)));
        [pd, j] = min(hs.q{i}.pixdim);
        hs.q{i}.Rvec = R0 / R * pd / hs.q{hs.iback}.pixdim(j);
        
        guidata(hs.fig, hs);
    end
end

img = hs.q{i}.nii.img;
% This is needed since vec is in image ref, at least for fsl
img(:,:,:,hs.q{i}.flip) = -img(:,:,:,hs.q{i}.flip);
if isfield(hs.q{i}, 'mask') % ignore modulation
    img = bsxfun(@times, img, hs.q{i}.mask);
end
if any(abs(diff(hs.q{hs.iback}.pixdim))>1e-4) % non isovoxel background
    pd = hs.q{hs.iback}.pixdim;
    pd = pd / min(pd);
    for j = 1:3, img(:,:,:,j) = img(:,:,:,j) / pd(j); end
end

if isfield(hs.q{i}, 'Rvec')
    img = reshape(img, [prod(d(1:3)) d(4)]);
    img = img * hs.q{i}.Rvec;
    img = reshape(img, d);
end

for ix = iaxis
    I = round(get(hs.ijk(ix), 'Value'));
    j = 1:3; j(ix) = [];
    if isfield(hs.q{i}, 'ivec')
        I = abs(hs.q{i}.ivec(:,:,:,ix) - I);
        [~, I] = min(I, [], ix);
        
        ii = {1:d(1) 1:d(2) 1:d(3)};
        ii{ix} = single(1);
        [ii{1}, ii{2}, ii{3}] = ndgrid(ii{:});
        ii{ix} = single(I);
        io = {':' ':' ':' ':'}; io{ix} = 1;
        im = img(io{:});
        for k = 1:2
            im(:,:,:,k) = interp3(img(:,:,:,j(k)), ii{[2 1 3]}, 'nearest');
        end
    
        ind = sub2ind(d(1:3), ii{:});
        X = hs.q{i}.ivec(:,:,:,j(1)); X = permute(X(ind), [j([2 1]) ix]);
        Y = hs.q{i}.ivec(:,:,:,j(2)); Y = permute(Y(ind), [j([2 1]) ix]);    
    else
        ii = {':' ':' ':'};
        ii{ix} = I;
        im = img(ii{:}, j);
        [Y, X] = ndgrid(1:d(j(2)), 1:d(j(1)));
    end
    
    im = permute(im, [j([2 1]) 4 ix]);
    im(im==0) = nan; % avoid dots in emf and eps
    X = X - im(:,:,1)/2;
    Y = Y - im(:,:,2)/2;
    set(p(i).hsI(ix), 'XData', X, 'YData', Y, 'UData', im(:,:,1), 'VData', im(:,:,2));
end

%% Bring cross and label to front
function crossFront(hs)
for i = 1:3
    txt = allchild(hs.ax(i));
    ind = strcmpi(get(txt, 'Type'), 'text');
    txt = txt(ind); % a number, two letters, plus junk text with matlab 2010b
    uistack([txt' hs.cross(i,:)], 'top');
end

%% Compute color map for LUT
function map = lut2map(lut, hs)
map = linspace(0,1,64)'*[1 1 1]; % gray
if     lut == 1, return; % gray
elseif lut == 2, map(:,2:3) = 0; % red
elseif lut == 3, map(:,[1 3]) = 0; % green
elseif lut == 4, map(:,1:2) = 0; % blue
elseif lut == 5, map(:,2) = 0; % violet
elseif lut == 6, map(:,3) = 0; % yellow
elseif lut == 7, map(:,1) = 0; % cyan
elseif any(lut == [8 19 26]), map(:,3) = 0; map(:,1) = 1; % red_yellow
elseif lut == 9, map(:,1) = 0; map(:,3) = map(end:-1:1,3); % blue_green
elseif lut == 10 % two-sided
    map = map(1:2:end,:); % half
    map_neg = map;
    map(:,3) = 0; map(:,1) = 1; % red_yellow
    map_neg(:,1) = 0; map_neg(:,3) = map_neg(end:-1:1,3); % blue_green
    map = [map_neg(end:-1:1,:); map];
elseif lut == 11, map(:,2:3) = 0; % vector lines
elseif lut == 12 % parula not in old matlab, otherwise this can be omitted
    map = [ 0.208 0.166 0.529
            0.212 0.190 0.578
            0.212 0.214 0.627
            0.208 0.239 0.677
            0.196 0.264 0.728
            0.171 0.292 0.779
            0.125 0.324 0.830
            0.059 0.360 0.868
            0.012 0.388 0.882
            0.006 0.409 0.883
            0.017 0.427 0.879
            0.033 0.443 0.872
            0.050 0.459 0.864
            0.063 0.474 0.855
            0.072 0.489 0.847
            0.078 0.504 0.838
            0.079 0.520 0.831
            0.075 0.538 0.826
            0.064 0.557 0.824
            0.049 0.577 0.823
            0.034 0.597 0.820
            0.026 0.614 0.814
            0.024 0.629 0.804
            0.023 0.642 0.791
            0.023 0.653 0.777
            0.027 0.664 0.761
            0.038 0.674 0.744
            0.059 0.684 0.725
            0.084 0.693 0.706
            0.113 0.702 0.686
            0.145 0.710 0.665
            0.180 0.718 0.642
            0.218 0.725 0.619
            0.259 0.732 0.595
            0.302 0.738 0.571
            0.348 0.742 0.547
            0.395 0.746 0.524
            0.442 0.748 0.503
            0.487 0.749 0.484
            0.530 0.749 0.466
            0.571 0.749 0.449
            0.610 0.747 0.434
            0.647 0.746 0.419
            0.683 0.743 0.404
            0.718 0.741 0.390
            0.752 0.738 0.377
            0.786 0.736 0.363
            0.819 0.733 0.350
            0.851 0.730 0.336
            0.882 0.727 0.322
            0.914 0.726 0.306
            0.945 0.726 0.289
            0.974 0.731 0.267
            0.994 0.745 0.240
            0.999 0.765 0.216
            0.996 0.786 0.197
            0.988 0.807 0.179
            0.979 0.827 0.163
            0.970 0.848 0.147
            0.963 0.871 0.131
            0.959 0.895 0.113
            0.960 0.922 0.095
            0.966 0.951 0.076
            0.976 0.983 0.054];
elseif lut < 26 % matlab LUT
    map = feval(hs.lutStr{lut}, 64);
elseif lut == 27 % phase3: red-yellow-green-yellow-red
    a = map(:,1);
    map(1:32,3) = 0; map(1:16,1) = 1; map(17:32,2) = 1;
    map(1:16,2) = a(1:4:64); map(17:32,1) = a(64:-4:1);
    map(33:64,:) = map(32:-1:1,:);
elseif lut == 28 % phase6: red-yellow-green/violet-blue-cyan
    a = map(:,1);
    map(1:32,3) = 0; map(1:16,1) = 1; map(17:32,2) = 1;
    map(1:16,2) = a(1:4:64); map(17:32,1) = a(64:-4:1);
    map(33:48,2) = 0; map(33:48,3) = 1; map(33:48,1) = a(64:-4:1);
    map(49:64,1) = 0; map(49:64,3) = 1; map(49:64,2) = a(1:4:64);
elseif lut == 29 % RGB
else % custom
    map = hs.lutStr{lut};
end

%% Preference dialog
function pref_dialog(pref)
pf = get(pref, 'UserData');
d = dialog('Name', 'Preferences', 'Visible', 'off');
pos = getpixelposition(d);
pos(3:4) = [396 332];
h.fig = ancestor(pref, 'figure');  

uicontrol(d, 'Style', 'text', 'Position', [8 306 300 22], ...
    'String', 'Background (template) image folder:', 'HorizontalAlignment', 'left');
h.openPath = uicontrol(d, 'Style', 'edit', 'String', pf.openPath, ...
    'Position', [8 288 350 22], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
    'TooltipString', 'nii_viewer will point to this folder when you "Open" image');
uicontrol('Parent', d, 'Position', [358 289 30 22], 'Tag', 'browse', ...
    'String', '...', 'Callback', @pref_dialog_cb);

h.rightOnLeft = uicontrol(d, 'Style', 'popup', 'BackgroundColor', 'w', ...
    'Position', [8 252 380 22], 'Value', pf.rightOnLeft+1, ...
    'String', {'Neurological orientation (left on left side)' ...
               'Radiological orientation (right on left side)'}, ...
    'TooltipString', 'Display convention also applies to future use');

uicontrol(d, 'Style', 'text', 'Position', [8 210 40 22], ...
    'String', 'Layout', 'HorizontalAlignment', 'left', ...
    'TooltipString', 'Layout for three views');

% iconsFolder = fullfile(matlabroot,'/toolbox/matlab/icons/');
% iconUrl = strrep(['file:/' iconsFolder 'matlabicon.gif'],'\','/');
% str = ['<html><img src="' iconUrl '"/></html>'];
h.layout = uicontrol(d, 'Style', 'popup', 'BackgroundColor', 'w', ...
    'Position', [50 214 338 22], 'Value', pf.layout, ...
    'String', {'one-row' 'two-row sag on right' 'two-row sag on left'}, ...
    'TooltipString', 'Layout for three views');

h.mouseOver = uicontrol(d, 'Style', 'checkbox', ...
    'Position', [8 182 380 22], 'Value', pf.mouseOver, ...
    'String', 'Show coordinates and intensity when mouse moves over image', ...
    'TooltipString', 'Also apply to future use');

uipanel(d, 'Units', 'Pixels', 'Position', [4 110 390 56], 'BorderType', 'etchedin', ...
    'BorderWidth', 2, 'Title', 'For "Save NIfTI as" if interpolation is applicable');
str = {'nearest' 'linear' 'cubic' 'spline'};
val = find(strcmp(str, pf.interp));
uicontrol('Parent', d, 'Style', 'text', 'Position', [8 116 140 22], ...
    'String', 'Interpolation method:', 'HorizontalAlignment', 'right');
h.interp = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [150 120 68 22], 'Value', val, 'BackgroundColor', 'w');

uicontrol('Parent', d, 'Style', 'text', 'Position', [230 116 90 22], ...
    'String', 'Missing value:', 'HorizontalAlignment', 'right');
h.extraV = uicontrol(d, 'Style', 'edit', 'String', num2str(pf.extraV), ...
    'Position', [324 120 60 22], 'BackgroundColor', 'w', ...
    'TooltipString', 'NaN or 0 is typical, but can be any number');

str = strtrim(cellstr(num2str([0 120 150 200 300 600 1200]')));
val = find(strcmp(str, pf.dpi));
uipanel(d, 'Units', 'Pixels', 'Position', [4 40 390 56], 'BorderType', 'etchedin', ...
    'BorderWidth', 2, 'Title', 'For "Save figure as" and "Copy figure"');
uicontrol('Parent', d, 'Style', 'text', 'Position', [8 46 90 22], ...
    'String', 'Resolution:', 'HorizontalAlignment', 'right');
h.dpi = uicontrol(d, 'Style', 'popup', 'String', str, ...
    'Position', [110 50 50 22], 'Value', val, 'BackgroundColor', 'w', ...
    'TooltipString', 'in DPI (0 means screen resolution)');

uicontrol('Parent', d, 'Position', [300 10 70 24], 'Tag', 'OK', ...
    'String', 'OK', 'Callback', {@pref_dialog_cb, pref});
uicontrol('Parent', d, 'Position',[200 10 70 24], ...
    'String', 'Cancel', 'Callback', 'delete(gcf)');

set(d, 'Position', pos, 'Visible', 'on');
guidata(d, h);

%% Preference dialog callback
function pref_dialog_cb(h, ~, pref)
hs = guidata(h);
if strcmp(get(h, 'Tag'), 'OK') % done
    pf = get(pref, 'UserData');
    fh = hs.fig;
    
    pf.rightOnLeft = get(hs.rightOnLeft, 'Value')==2;
    hLR = findobj(fh, 'Type', 'uimenu', 'Label', 'Right on left side');
    if strcmp(get(hLR, 'Checked'), 'off') == pf.rightOnLeft
        nii_viewer_cb(hLR, [], 'flipLR', fh);
    end
    
    pf.mouseOver = get(hs.mouseOver, 'Value');
    if pf.mouseOver
        set(fh, 'WindowButtonMotionFcn', {@nii_viewer_cb 'mousemove' fh});
    else
        set(fh, 'WindowButtonMotionFcn', '');        
    end
    
    pf.layout = get(hs.layout, 'Value'); % 0 or 1
    
    i = get(hs.interp, 'Value');
    str = get(hs.interp, 'String');
    pf.interp = str{i};
    
    pf.extraV = str2double(get(hs.extraV, 'String'));
    pf.openPath = get(hs.openPath, 'String');

    i = get(hs.dpi, 'Value');
    str = get(hs.dpi, 'String');
    pf.dpi = str{i}; 
        
    set(pref, 'UserData', pf);
    setpref('nii_viewer_para', fieldnames(pf), struct2cell(pf));
    delete(get(h, 'Parent'));
elseif strcmp(get(h, 'Tag'), 'browse') % set openPath
    pth = uigetdir(pwd, 'Select folder for background image');
    if ~ischar(pth), return; end
    set(hs.openPath, 'String', pth);
end

%% Simple version of interp3
function V = interp3a(V, I, method)
% V = interp3a(V, I, 'linear');
% This is similar to interp3 from Matlab, but takes care of the Matlab version
% issue, and the input is simplified for coordinate. The coordinate input are in
% this way: I(1,:), I(2,:) and I(3,:) are for x, y and z. 
persistent v;
if isempty(v)
    try 
        griddedInterpolant({1:3, 1:3, 1:3}, ones(3,3,3), 'nearest', 'none');
        v = 2014;
    catch
        try
            griddedInterpolant({1:3, 1:3, 1:3}, ones(3,3,3), 'nearest');
            v = 2013;
        catch
            v = 2011;
        end
    end
end
if v > 2011
    d = size(V); d(numel(d)+1:3) = 1;
    if  v > 2013
        F = griddedInterpolant({1:d(1), 1:d(2), 1:d(3)}, V, method, 'none');
    else
        F = griddedInterpolant({1:d(1), 1:d(2), 1:d(3)}, V, method);
    end
    V = F(I(1,:), I(2,:), I(3,:)); % interpolate
else % earlier matlab
    V = interp3(V, I(2,:), I(1,:), I(3,:), method, nan);
end

%% 2D/3D smooth wrapper: no input check for 2D
function out = smooth23(in, method, n, varargin)
% out = smooth23(in, method, n, varargin)
% This works the same as smooth3 from Matlab, except it also works if the input
% is 2D.
if size(in,3)>1, out = smooth3(in, method, n, varargin{:}); return; end
if nargin<3 || isempty(n), n = 3; end
if numel(n)==1, n = [n n]; end
k = floor(n/2);
if k<1, out = in; return; end
n = k * 2 + 1; % odd number
if strcmp(method, 'box')
    kernal = ones(n) / n(1)/n(2);
elseif strcmp(method, 'gaussian')
    if nargin<4 || isempty(varargin{1}), sd = 0.65; 
    else, sd = varargin{1}; 
    end
    [x, y] = ndgrid(-k(1):k(1), -k(2):k(2));
    kernal = exp(-(x.*x  +  y.*y) / (2*sd*sd));
    kernal = kernal / sum(kernal(:));
else
    fprintf(2, 'Invalid smooth method: %s\n', method);
    out = in; return;
end
in = [repmat(in(1,:), [k(1) 1]); in; repmat(in(end,:), [k(1) 1])];
in = [repmat(in(:,1), [1 k(2)])  in  repmat(in(:,end), [1 k(2)])];
out = conv2(in, kernal, 'valid');

%% Show ijk/xyz and value
function xyz = set_xyz(hs, I)
if nargin<2
    for i=3:-1:1, I(i) = get(hs.ijk(i), 'Value'); end
end

p = hs.scroll.UserData;
I = round(I);
xyz = round(hs.q{hs.iback}.R * [I-1 1]');
xyz = xyz(1:3);
str = sprintf('(%g,%g,%g)=(%i,%i,%i): ', I, xyz);

for i = 1:numel(p) % show top one first
    if p(i).show == 0, continue; end
    t = round(p(i).volume);
    if isfield(hs.q{i}, 'R0') % need interpolation
        if isfield(hs.q{i}, 'warp')
            warp = hs.q{i}.warp(I(1), I(2), I(3), :);
            I0 = hs.q{i}.Ri * (hs.q{i}.R0 * [I-1 1]' + [warp(:); 0]);
        else
            I0 = hs.q{i}.Ri * hs.q{i}.R0 * [I-1 1]'; % overlay ijk
        end
        I0 = round(I0(1:3)+1);
    else, I0 = I;
    end
    try
        val = hs.q{i}.nii.img(I0(1), I0(2), I0(3), t, :);
        if isfield(hs.q{i}, 'scl_slope')
            val = val * hs.q{i}.scl_slope + hs.q{i}.scl_inter;
        end
    catch
        val = nan; % out of range
    end
    
    if isfield(hs.q{i}, 'labels')
        try 
            labl = hs.q{i}.labels{val};
            str = sprintf('%s %s', str, labl);
            continue; % skip numeric val assignment
        catch
        end
    end
    
    fmtstr = '%.5g ';
    if numel(val)>1
        fmtstr = repmat(fmtstr, 1, numel(val));
        fmtstr = ['[' fmtstr]; fmtstr(end) = ']'; %#ok
    end
    str = sprintf(['%s ' fmtstr], str, val);
end
hs.value.String = str;

%% nii essentials
function s = nii_essential(hdr)
% info = nii_essential(hdr);
% Decode important NIfTI hdr into struct info, which is human readable.
if isfield(hdr, 'nii') % q input by nii_viewer
    s.FileName = hdr.nii.hdr.file_name;
    if isfield(hdr, 'mask_info'), s.MaskedBy = hdr.mask_info; end
    if isfield(hdr, 'modulation_info'), s.ModulatdBy = hdr.modulation_info; end
    if isfield(hdr, 'alignMtx'), s.AlignMatrix = hdr.alignMtx; end
    hdr = hdr.nii.hdr;
else
    s.FileName = hdr.file_name;
end
switch hdr.intent_code
    case 2, s.intent = 'Correlation'; s.DoF = hdr.intent_p1;
    case 3, s.intent = 'T-test';      s.DoF = hdr.intent_p1;
    case 4, s.intent = 'F-test';      s.DoF = [hdr.intent_p1 hdr.intent_p2];
    case 5, s.intent = 'Z-score';
    case 6, s.intent = 'Chi squared'; s.DoF = hdr.intent_p1;
        
    % several non-statistical intent_code
    case 1002, s.intent = 'Label'; % e.g. AAL labels
    case 2003, s.intent = 'RGB'; % triplet in the 5th dimension
    case 2004, s.intent = 'RGBA'; % quadruplet in the 5th dimension
end
switch hdr.datatype
    case 0
    case 1,    s.DataType = 'logical';
    case 2,    s.DataType = 'uint8';
    case 4,    s.DataType = 'int16';
    case 8,    s.DataType = 'int32';
    case 16,   s.DataType = 'single';
    case 32,   s.DataType = 'single complex';
    case 64,   s.DataType = 'double';
    case 128,  s.DataType = 'uint8 RGB';
    case 256,  s.DataType = 'int8';
    case 511,  s.DataType = 'single RGB';
    case 512,  s.DataType = 'uint16';
    case 768,  s.DataType = 'uint32';
    case 1024, s.DataType = 'int64';
    case 1280, s.DataType = 'uint64';
    case 1792, s.DataType = 'double complex';
    case 2304, s.DataType = 'uint8 RGBA';
    otherwise, s.DataType = 'unknown';
end
s.Dimension = hdr.dim(2:hdr.dim(1)+1);
switch bitand(hdr.xyzt_units, 7)
    case 1, s.VoxelUnit = 'meters';
    case 2, s.VoxelUnit = 'millimeters';
    case 3, s.VoxelUnit = 'micrometers';
end 
s.VoxelSize = hdr.pixdim(2:4);
switch bitand(hdr.xyzt_units, 56)
    case 8,  s.TemporalUnit = 'seconds';
    case 16, s.TemporalUnit = 'milliseconds';
    case 24, s.TemporalUnit = 'microseconds';
    case 32, s.TemporalUnit = 'Hertz';
    case 40, s.TemporalUnit = 'ppm';
    case 48, s.TemporalUnit = 'radians per second';
end 
if isfield(s, 'TemporalUnit') && strfind(s.TemporalUnit, 'seconds')
    s.TR = hdr.pixdim(5);
end
if hdr.dim_info>0
    s.FreqPhaseSliceDim = bitand(hdr.dim_info, [3 12 48]) ./ [1 4 16];
    a = bitand(hdr.dim_info, 192) / 64;
    if a>0 && s.FreqPhaseSliceDim(2)>0
        ax = 'xyz'; % ijk to be accurate
        pm = ''; if a == 2, pm = '-'; end 
        s.PhaseEncodingDirection = [pm ax(s.FreqPhaseSliceDim(2))]; 
    end
end

switch hdr.slice_code
    case 0
    case 1, s.SliceOrder = 'Sequential Increasing';
    case 2,	s.SliceOrder = 'Sequential Decreasing';
    case 3,	s.SliceOrder = 'Alternating Increasing 1';
    case 4,	s.SliceOrder = 'Alternating Decreasing 1';
    case 5,	s.SliceOrder = 'Alternating Increasing 2';
    case 6,	s.SliceOrder = 'Alternating Decreasing 2';
    otherwise, s.SliceOrder = 'Multiband?';
end
if ~isempty(hdr.descrip), s.Notes = hdr.descrip; end
str = formcode2str(hdr.qform_code);
if ~isempty(str), s.qform = str; end
if hdr.qform_code>0, s.qform_mat = quat2R(hdr); end
str = formcode2str(hdr.sform_code);
if ~isempty(str), s.sform = str; end
if hdr.sform_code>0
    s.sform_mat = [hdr.srow_x; hdr.srow_y; hdr.srow_z; 0 0 0 1];
end

%% decode NIfTI form_code
function str = formcode2str(code)
switch code
    case 0, str = '';
    case 1, str = 'Scanner Anat';
    case 2, str = 'Aligned Anat';
    case 3, str = 'Talairach';
    case 4, str = 'mni_152';
    otherwise, str = 'Unknown';
end

%% Get a mask based on image intensity, but with inside brain filled
function r = img2mask(img)
mn = mean(img(img(:)>0));
r = smooth23(img, 'box', 5) > mn/8; % smooth, binarize
if sum(r(:))==0, return; end

try
    C = contourc(double(r), [1 1]);
    i = 1; c = {};
    while size(C,2)>2 % split C into contours
        k = C(2,1) + 1;
        c{i} = C(:, 2:k); C(:,1:k) = []; %#ok
        i = i+1;
    end
    
    nc = numel(c);
    rg = nan(nc, 4); % minX minY maxX maxY
    for i = 1:nc
        rg(i,:) = [min(c{i},[],2)' max(c{i},[],2)'];
    end
    ind = false(nc,1);
    foo = min(rg(:,1)); ind = ind | foo==rg(:,1);
    foo = min(rg(:,2)); ind = ind | foo==rg(:,2);
    foo = max(rg(:,3)); ind = ind | foo==rg(:,3);
    foo = max(rg(:,4)); ind = ind | foo==rg(:,4);
    c = c(ind); % outmost contour(s) 
    len = cellfun(@(x) size(x,2), c);
    [~, ind] = sort(len, 'descend');
    c = c(ind);
    C = c{1};
    if isequal(C(:,1), C(:,end)), c(2:end) = []; end % only 1st if closed
    nc = numel(c);
    for i = nc:-1:2 % remove closed contours except one with max len
        if isequal(c{i}(:,1), c{i}(:,end)), c(i) = []; end
    end
    nc = numel(c);
    while nc>1 % +1 contours, put all into 1st
        d2 = nan(nc-1, 2); % distance^2 from C(:,end) to other start/endpoint
        for i = 2:nc
            d2(i-1,:) = sum((C(:,end)*[1 1] - c{i}(:,[1 end])).^2);
        end
        [i, j] = find(d2 == min(d2(:)));
        i = i + 1; % start with 2nd
        if j == 1, C = [C c{i}]; %#ok C(:,1) connect to c{i}(:,1)
        else C = [C c{i}(:,end:-1:1)]; %#ok C(:,end) to c{i}(:,end)
        end
        c(i) = []; nc = nc-1;
    end
    if ~isequal(C(:,1), C(:,end)), C(:,end+1) = C(:,1); end % close the contour
    x = C(1, :);
    y = C(2, :);
    [m, n] = size(r);

    % following is the method in Octave poly2mask
    xe = [x(1:numel(x)-1); x(1, 2:numel(x))]; % edge x
    ye = [y(1:numel(y)-1); y(1, 2:numel(y))]; % edge y
    ind = ye(1,:) == ye(2,:);
    xe(:,ind) = []; ye(:, ind) = []; % reomve horizontal edges
    minye = min(ye);
    maxye = max(ye);
    t = (ye == [minye; minye]);
    exminy = xe(:); exminy = exminy(t);
    exmaxy = xe(:); exmaxy = exmaxy(~t);
    maxye = maxye';
    minye = minye';
    m_inv = (exmaxy - exminy) ./ (maxye - minye);
    ge = [maxye minye exmaxy m_inv];
    ge = sortrows(ge, [1 3]);
    ge = [-Inf -Inf -Inf -Inf; ge];

    gei = size(ge, 1);
    sl = ge(gei, 1);
    ae = []; % active edge
    while (sl == ge(gei, 1))
        ae = [ge(gei, 2:4); ae]; %#ok
        gei = gei - 1;
    end

    miny = min(y);
    if miny < 1, miny = 1; end

    while (sl >= miny)
        if (sl <= m) % check vert clipping
            ie = round(reshape(ae(:, 2), 2, size(ae, 1)/2));
            ie(1, :) = ie(1, :) + (ie(1, :) ~= ie(2, :));
            ie(1, (ie(1, :) < 1)) = 1;
            ie(2, (ie(2, :) > n)) = n;
            ie = ie(:, (ie(1, :) <= n));
            ie = ie(:, (ie(2, :) >= 1));
            for i = 1:size(ie,2)
                r(sl, ie(1, i):ie(2, i)) = true;
            end
        end

        sl = sl - 1;
        ae = ae((ae(:, 1) ~= sl), :);
        ae(:, 2) = ae(:, 2) - ae(:, 3);

        while (sl == ge(gei, 1))
            ae = [ae; ge(gei, 2:4)]; %#ok
            gei = gei - 1;
        end

        if size(ae,1) > 0
            ae = sortrows(ae, 2);
        end
    end
catch %me, fprintf(2, '%s\n', me.message); assignin('base', 'me', me);
end

%% update colorbar label
function set_colorbar(hs)
if strcmpi(get(hs.colorbar, 'Visible'), 'off'), return; end
p = hs.scroll.UserData(hs.files.getSelectedIndex+1);
if p.lut == 11, map = lut2map(hs.lut.UserData, hs);
else, map = lut2map(p.lut, hs);
end
rg = sort([p.lb p.ub]);
if any(p.lut == 26:28)
    labls = [0 180 360];
elseif p.lut~=10
    if rg(2)<0, rg = rg([2 1]); end
    mn = str2double(num2str(mean(rg), '%.4g'));
    labls = [rg(1) mn rg(2)];
else
    rg = sort(abs(rg));
    labls = {num2str(-rg(2)) num2str(rg(1),'+/-%g') num2str(rg(2))};
end
% colormap in earlier matlab version changes values in colorbar.
% So we have to set those values each time.
% set(hs.colorbar, 'YTickLabel', labls); % new matlab
colormap(hs.ax(4), map);
set(get(hs.colorbar, 'Children'), 'YData', [0 1]); % Trick for old matlab
set(hs.colorbar, 'YTickLabel', labls, 'YTick', [0 0.5 1], 'Ylim', [0 1]);

%% return screen size in pixels
function res = screen_pixels(id)
res = get(0, 'MonitorPositions');
if size(res,1)<2, res = res(1, 3:4); return; end % single/duplicate monitor
if nargin, res = res(id,3:4); return; end
res = sortrows(res);
res = res(end,1:2) + res(end,3:4) - res(1,1:2);

%% add mask or modulation
function addMask(h, ~)
hs = guidata(h);
i = hs.files.getSelectedIndex+1;
pName = fileparts(hs.scroll.UserData(i).fname);
[fname, pName] = uigetfile([pName '/*.nii;*.hdr;*.nii.gz;*.hdr.gz'], ...
    'Select mask NIfTI');
if ~ischar(fname), return; end
fname = fullfile(pName, fname);

nii = nii_tool('load', fname);
hdr = hs.q{i}.nii.hdr;
codes = [hdr.sform_code hdr.qform_code];
[R, frm] = nii_xform_mat(nii.hdr, codes);
if ~any(frm == codes)
    str = ['There is no matching coordinate system between the selected ' ...
        'image and the mask image. Do you want to apply the mask anyway?'];
    btn = questdlg(str, 'Apply mask?', 'Cancel', 'Apply', 'Cancel');
    if isempty(btn) || strcmp(btn, 'Cancel'), return; end
    R0 = nii_xform_mat(hdr, codes(1));
else
    R0 = nii_xform_mat(hdr, frm); % may be the same as hs.q{i}.R
end
R0 = reorient(R0, hdr.dim(2:4)); % do this since img was done when loaded

% if isfield(hs.q{i}, 'alignMtx'), R = R0 / hs.q{i}.R * R; end % inverse
% this wont work if lines is background & Mprage is overlay
if all(isfield(hs.q{i}, {'R0' 'alignMtx'})) % target as mask
    R1 = reorient(R, nii.hdr.dim(2:4));
    if all(abs(R1(:)-hs.q{i}.R0(:))<1e-4), R0 = hs.q{i}.R; end % not 100% safe
end

d = single(size(hs.q{i}.nii.img)); % dim for reoriented img
d(numel(d)+1:3) = 1; d = d(1:3);

I = ones([d 4], 'single');
[I(:,:,:,1), I(:,:,:,2), I(:,:,:,3)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
I = permute(I, [4 1 2 3]);
I = reshape(I, [4 prod(d)]); % ijk grids of target img
I = inv(R) * R0 * I + 1; %#ok ijk+1 for mask
I = round(I * 100) / 100;

im = single(nii.img(:,:,:,1)); % first mask volume
slope = nii.hdr.scl_slope;
if slope==0, slope = 1; end
im = im * slope + nii.hdr.scl_inter;
im = interp3a(im, I, 'nearest');
im1 = im(~isnan(im)); % for threshold computation
im = reshape(im, d);

if strcmp(get(h, 'Label'), 'Apply mask') % binary mask
    if numel(unique(im1))<3
        thre = min(im1);
    else
        a = get_range(nii);
        str = sprintf('Threshold for non-binary mask (%.3g to %.4g)', ...
            min(im1), max(im1));
        a = inputdlg(str, 'Input mask threshold', 1, {num2str(a(1), '%.3g')});
        if isempty(a), return; end
        thre = str2double(a{1});
        fname = [fname ' (threshold = ' a{1} ')']; % info only
    end
    hs.q{i}.mask = ones(size(im), 'single');
    hs.q{i}.mask(abs(im)<=thre) = nan;
    hs.q{i}.mask_info = fname;
    noteStr = '(masked)';
else % modulation
    mi = min(im1); ma = max(im1);
    if mi<0 || ma>1
        str = {sprintf('Lower bound to clip to 0 (image min = %.2g)', mi) ...
               sprintf('Upper bound to clip to 1 (image max = %.2g)', ma)};
        def = strtrim(cellstr(num2str([mi;ma], '%.2g'))');
        a = inputdlg(str, 'Input modulation range', 1, def);
        if isempty(a), return; end
        mi = str2double(a{1}); ma = str2double(a{2});
        fname = [fname ' (range [' a{1} ' ' a{2} '])'];
    end
    im(im<mi) = mi; im(im>ma) = ma;
    hs.q{i}.modulation = (im-mi) / (ma-mi);
    hs.q{i}.modulation_info = fname;
    noteStr = '(modulated)';
end
guidata(hs.fig, hs);
set_cdata(hs);

str = hs.files.getModel.get(i-1);
n = numel(noteStr);
if numel(str)<n || ~strcmp(str(end+(-n+1:0)), noteStr)
    hs.files.getModel.set(i-1, [str noteStr]);
end

%% update crosshair: ix correspond to one of the three spinners, not views
function set_cross(hs, ix)
h = hs.cross;
for i = ix
    c = get(hs.ijk(i), 'Value');
    g = hs.gap(i);
    if i == 1 % I changed
        set([h(2,3:4) h(3,3:4)], 'XData', [c c]);
        set([h(2,1) h(3,1)], 'XData', [0 c-g]);
        set([h(2,2) h(3,2)], 'XData', [c+g hs.dim(1)+1]);
    elseif i == 2 % J
        set(h(1,3:4), 'XData', [c c]);
        set(h(1,1), 'XData', [0 c-g]);
        set(h(1,2), 'XData', [c+g hs.dim(2)+1]);
        set(h(3,1:2), 'YData', [c c]);
        set(h(3,3), 'YData', [0 c-g]);
        set(h(3,4), 'YData', [c+g hs.dim(2)+1]);
    else % K
        set([h(1,1:2) h(2,1:2)], 'YData', [c c]);
        set([h(1,3) h(2,3)], 'YData', [0 c-g]);
        set([h(1,4) h(2,4)], 'YData', [c+g hs.dim(3)+1]);
    end
end

%% Switch display parameters to selected file
function set_file(~, evt, fh)
if ~isempty(evt) && evt.getValueIsAdjusting, return; end
hs = guidata(fh);
i = hs.files.getSelectedIndex+1;
if i<1 || i>numel(hs.q), return; end
p = hs.scroll.UserData(i);
nam = {'lb' 'ub' 'alpha' 'volume' 'lut' 'smooth' 'interp' };
cb = cell(4,1);
for j = 1:4 % avoid firing spinner callback
    cb{j} = hs.(nam{j}).StateChangedCallback;
    hs.(nam{j}).StateChangedCallback = '';
end

for j = 1:numel(nam)
    set(hs.(nam{j}), 'Value', p.(nam{j}));
end
set(hs.lb.Model, 'StepSize', p.lb_step);
set(hs.ub.Model, 'StepSize', p.ub_step);
nVol = size(hs.q{i}.nii.img, 4);
set(hs.volume, 'Enable', nVol>1, 'ToolTipText', ['Volume number, 1:' num2str(nVol)]);
set(hs.volume.Model, 'Maximum', nVol);

for j = 1:4 % restore spinner callback
    hs.(nam{j}).StateChangedCallback = cb{j};
end
set_colorbar(hs);

off_on = {'off' 'on'};
set(hs.interp, 'Enable', off_on{isfield(hs.q{i}, 'R0')+1});
set(hs.overlay(1:4), 'Enable', off_on{(numel(hs.q)>1)+1}); % stack & Close overlays
set(hs.overlay(5), 'Enable', off_on{(i~=hs.iback)+1}); % Close overlay

%% Duplicate image handles, inlcuding ButtonDownFcn for new matlab
function h = copyimg(hs)
h = hs.hsI;
for i = 1:3, h(i) = copyobj(hs.hsI(i), hs.ax(i)); end
cb = get(hs.hsI(1), 'ButtonDownFcn');
set(h, 'Visible', 'on', 'ButtonDownFcn', cb); % cb needed for 2014b+
crossFront(hs);

%% Save selected nii as another file
function save_nii_as(h, ~)
hs = guidata(h);
c = get(h, 'Label');
i = hs.files.getSelectedIndex+1;
nam = hs.scroll.UserData(i).fname;
pName = fileparts(nam);
if isempty(pName), pName = pwd; end
try nii = nii_tool('load', nam); % re-load to be safe
catch % restore reoriented img
    nii = hs.q{i}.nii;
    if ~exist('flip', 'builtin'), eval('flip=@flipdim;'); end
    for k = 1:3, if hs.q{i}.flip(k), nii.img = flip(nii.img, k); end; end
    nii.img = permute(nii.img, [hs.q{i}.perm 4:8]); % all vol in dim(4)
    if nii.hdr.datatype == 4 % leave others as it is or single
        nii.img = int16(nii.img);
    elseif nii.hdr.datatype == 512
        nii.img = uint16(nii.img);
    elseif any(nii.hdr.datatype == [2 128 2304])
        nii.img = uint8(nii.img);
    end
end

if ~isempty(strfind(c, 'a copy')) %#ok<*STREMP> % a copy
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input name for FSL RGB file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname);
elseif ~isempty(strfind(c, 'dim 4')) % fsl RGB
    if any(size(nii.img,8) == 3:4)
        nii.img = permute(nii.img, [1:3 8 4:7]);
    elseif ~any(nii.hdr.dim(5) == 3:4)
        errordlg('Selected image is not RGB data.'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input name for FSL RGB file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname);
elseif ~isempty(strfind(c, 'dim 3')) % old mricron RGB
    if any(nii.hdr.dim(5) == 3:4)
        nii.img = permute(nii.img, [1:3 5:7 4]);
    elseif ~any(size(nii.img,8) == 3:4)
        errordlg('Selected image is not RGB data'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input name for old mricrom styte file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    old = nii_tool('RGBStyle', 'mricron');
    nii_tool('save', nii, fname);
    nii_tool('RGBStyle', old);
elseif ~isempty(strfind(c, 'AFNI')) % NIfTI RGB
    if any(nii.hdr.dim(5) == 3:4)
        nii.img = permute(nii.img, [1:3 5:8 4]);
    elseif ~any(size(nii.img,8) == 3:4)
        errordlg('Selected image is not RGB data'); return;
    end
    nii.img = abs(nii.img);
    [fname, pName] = uiputfile([pName '/*.nii'], ...
        'Input name for NIfTI standard RGB file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    old = nii_tool('RGBStyle', 'afni');
    nii_tool('save', nii, fname);
    nii_tool('RGBStyle', old);
elseif ~isempty(strfind(c, '3D')) % SPM 3D
    if nii.hdr.dim(5)<2
        errordlg('Selected image is not multi-volume data'); return;
    end
    [fname, pName] = uiputfile([pName '/*.nii'], 'Input base name for SPM 3D file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    nii_tool('save', nii, fname, 1); % force 3D
elseif ~isempty(strfind(c, 'new resolution'))
    str = 'Resolution for three dimension in mm';
    a = inputdlg(str, 'Input spatial resolution', 1, {'3 3 3'});
    if isempty(a), return; end
    res = sscanf(a{1}, '%g %g %g');
    if numel(res) ~= 3
        errordlg('Invalid spatial resolution');
        return;
    end
    if isequal(res, hs.q{i}.nii.hdr.pixdim(2:4))
        warndlg('The input resolution is the same as current one');
        return;
    end
    [fname, pName] = uiputfile([pName '/*.nii;nii.gz'], ...
        'Input result name for the new resolution file');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    pf = get(hs.pref, 'UserData');
    nii_xform(nii, res, fname, pf.interp, pf.extraV)
elseif ~isempty(strfind(c, 'matching background'))
    if i == hs.iback
        errordlg('You selected background image');
        return;
    end
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input result file name');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    pf = get(hs.pref, 'UserData');
    nii_xform(nii, hs.scroll.UserData(hs.iback).fname, fname, pf.interp, pf.extraV)
elseif ~isempty(strfind(c, 'aligned template'))
    [temp, pName] = uigetfile([pName '/*.nii;*.nii.gz'], ...
        'Select the aligned template file');
    if ~ischar(temp), return; end
    temp = fullfile(pName, temp);
    [mtx, pName] = uigetfile([pName '/*.mat'], ['Select the text ' ...
        'matrix file which aligns the nii to the template']);
    if ~ischar(mtx), return; end
    mtx = fullfile(pName, mtx);
    [fname, pName] = uiputfile([pName '/*.nii;*.nii.gz'], ...
        'Input result file name');
    if ~ischar(fname), return; end
    fname = fullfile(pName, fname);
    pf = get(hs.pref, 'UserData');
    nii_xform(nii, {temp mtx}, fname, pf.interp, pf.extraV)
else
    errordlg(sprintf('%s not implemented yet.', c));
end

%% Return 3-layer RGB, called by set_cdata
function [im, alfa] = lut2img(im, lut, rg, lutStr)
if any(lut == 26:28)
    im = im(:,:,2) .* single(im(:,:,1)>min(abs(rg))); % mag as mask
    rg = [0 1]; % disable later scaling
elseif isnumeric(lutStr) % custom LUT
    rg = [0 255]; % this seems right for aal and brodmann
end
if rg(2)<0 % asking for negative data
    rg = -rg([2 1]);
    if lut~=10, im = -im; end
end
if lut == 10 % two-sided, store negative value
    rg = sort(abs(rg));
    im_neg = -single(im) .* (im<0);
    im_neg = (im_neg-rg(1)) / (rg(2)-rg(1));
    im_neg(im_neg>1) = 1; im_neg(im_neg<0) = 0;
    alfa = im_neg; % add positive part later
    im_neg = repmat(im_neg, [1 1 3]); % gray now
else
    alfa = single(0);
end

if lut ~= 29 % not forced into RGB
    im = (im-rg(1)) / (rg(2)-rg(1));
    im(im>1) = 1; im(im<0) = 0;
    alfa = im + alfa;
    im = repmat(im, [1 1 3]); % gray now
end

switch lut
    case 1 % gray do nothing
    case 2, im(:,:,2:3) = 0; % red
    case 3, im(:,:,[1 3]) = 0; % green
    case 4, im(:,:,1:2) = 0; % blue
    case 5, im(:,:,2) = 0; % violet
    case 6, im(:,:,3) = 0; % yellow
    case 7, im(:,:,1) = 0; % cyan
    case {8 19} % red_yellow, autumn
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        im(:,:,3) = 0;
    case 9 % blue_green
        im(:,:,1) = 0;
        a = im(:,:,3); a(a==0) = 1; a = 1 - a; im(:,:,3) = a;
    case 10 % two-sided: combine red_yellow & blue_green
        im(:,:,3) = 0;
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        im_neg(:,:,1) = 0;
        a = im_neg(:,:,3); a(a==0) = 1; a = 1 - a; im_neg(:,:,3) = a;
        im = im + im_neg;
    case 15 % hot
        a = im(:,:,1); a = a/0.375; a(a>1) = 1; im(:,:,1) = a;
        a = im(:,:,2); a = a/0.375-1;
        a(a<0) = 0; a(a>1) = 1; im(:,:,2) = a;
        a = im(:,:,3); a = a*4-3; a(a<0) = 0; im(:,:,3) = a;
    case 16 % cool
        a = im(:,:,2); a(a==0) = 1; a = 1 - a; im(:,:,2) = a;
        a = im(:,:,3); a(a>0) = 1; im(:,:,3) = a;
    case 17 % spring
        a = im(:,:,1); a(a>0) = 1; im(:,:,1) = a;
        a = im(:,:,3); a(a==0) = 1; a = 1 - a; im(:,:,3) = a;
    case 18 % summer
        a = im(:,:,2); a(a==0) = -1; a = a/2+0.5; im(:,:,2) = a;
        a = im(:,:,3); a(a>0) = 0.4; im(:,:,3) = a;
    case 20 % winter
        im(:,:,1) = 0;
        a = im(:,:,3); a(a==0) = 2; a = 1-a/2; im(:,:,3) = a;
    case 22 % copper
        a = im(:,:,1); a = a*1.25; a(a>1) = 1; im(:,:,1) = a;
        im(:,:,2) = im(:,:,2) * 0.7812;
        im(:,:,3) = im(:,:,3) * 0.5;
    case 26 % phase, like red_yellow
        im(:,:,1) = 1; im(:,:,3) = 0;
    case 27 % phase3, red-yellow-green-yellow-red
        a = im(:,:,1);
        b1 = a<=0.25;
        b2 = a>0.25 & a<=0.5;
        b3 = a>0.5 & a<=0.75;
        b4 = a>0.75;
        a(b1 | b4) = 1;
        a(b2) = (0.5-a(b2))*4;
        a(b3) = (a(b3)-0.5)*4;
        im(:,:,1) = a;
        
        a = im(:,:,2);
        a(b2 | b3) = 1;
        a(b1) = a(b1)*4;
        a(b4) = (1-a(b4))*4;
        im(:,:,2) = a;
        
        im(:,:,3) = 0;
    case 28 % phase6, red-yellow-green/violet-blue-cyan
        a = im(:,:,1);
        b1 = a<=0.25;
        b2 = a>0.25 & a<=0.5;
        b3 = a>0.5 & a<=0.75;
        b4 = a>0.75;
        a(b2) =  (0.5-a(b2))*4; a(b1) = 1;
        a(b3) = (0.75-a(b3))*4; a(b4) = 0;
        im(:,:,1) = a;
        
        a = im(:,:,2);
        a(b1) = a(b1)*4; a(b2) = 1;
        a(b3) = 0; a(b4) = (a(b4)-0.75)*4;
        im(:,:,2) = a;
        
        a = im(:,:,3);
        a(b1 | b2) = 0;
        a(b3 | b4) = 1;
        im(:,:,3) = a;
    case 29 % disp non-NIfTI RGB as RGB
        if min(im(:)) < 0 % use abs if there is negative value
            rg(1) = 0;
            im = abs(im);
        end
        if max(im(:)) > 1
            im = (im-rg(1)) / (rg(2)-rg(1));
        end
        im(im>1) = 1; im(im<0) = 0;
        alfa = sum(im,3)/3;
    otherwise % parula(12), jet(13), hsv(14), bone(21), pink(23), custom
        if lut <= 25
            try
                map = feval(lutStr, 256);
            catch me
                if lut == 12 % parula is not in old matlab
                    map = lut2map(lut);
                else
                    rethrow(me);
                end
            end
            a = round(im(:,:,1) * size(map,1));
        else % custom
            map = lutStr; % LUT for custom actually
            a = round(im(:,:,1) * (size(map,1)-1)) + 1; % 1st is bkgrnd
        end
        
        for j = 1:size(a,1)
            for k = 1:size(a,2)
                if a(j,k)>0, im(j,k,:) = map(a(j,k),:); end
            end
        end
end

%% Return binary sphere ROI from xyz and r (mm)
function b = xyzr2roi(c, r, hdr)
% ROI_img = xyzr2roi(center, radius, hdr)
% Return an ROI img based on the dim info in NIfTI hdr. The center and radius
% are in unit of mm. 
d = single(hdr.dim(2:4));
I = ones([d 4], 'single');
[I(:,:,:,1), I(:,:,:,2), I(:,:,:,3)] = ndgrid(0:d(1)-1, 0:d(2)-1, 0:d(3)-1);
I = permute(I, [4 1 2 3]);
I = reshape(I, [4 prod(d)]); % ijk in 4 by nVox
R = nii_xform_mat(hdr);
I = R * I; % xyz in 4 by nVox

I = bsxfun(@minus, I(1:3,:), c(:)); % dist in x y z direction from center
I = sum(I .* I); % dist to center squared, 1 by nVox

b = I <= r*r; % within sphere
b = reshape(b, d);

%% Return center of gravity of an image
function c = img_cog(img)
% center_ijk = img_cog(img)
% Return the index of center of gravity in img (must be 3D).
img(isnan(img)) = 0;
img = double(abs(img));
gs = sum(img(:));
c = ones(3,1);
for i = 1:3
    if size(img,i)==1, continue; end
    a = shiftdim(img, i-1);
    a = sum(sum(a,3),2);
    c(i) = (1:size(img,i)) * a / gs;
end

%% set up disp parameter for new nifti in p(1)
function p = dispPara(p, hdr)
p(1).fname = hdr.file_name;
p(1).show = true; % img on
if any(hdr.datatype == [32 1792]) % complex
    p(1).lut = 26; % phase
    p(1).lb = str2double(sprintf('%.2g', p(1).ub/2));
elseif hdr.intent_code == 1002 % Label
    p(1).lut = 24; % prism
elseif hdr.intent_code > 0 % some stats
    if p(1).lb < 0
        p(1).lb = str2double(sprintf('%.2g', p(1).ub/2));
        p(1).lut = 10; % two-sided
    else
        a = setdiff(7:8, [p.lut]); % red-yellow & blue-green
        if isempty(a), a = 7; end % red-yellow
        p(1).lut = a(1);
    end
elseif numel(p) < 2
    p(1).lut = 1; % gray
else
    a = setdiff(2:7, [p.lut]); % use smallest unused mono-color lut 
    if isempty(a), a = 2; end % red
    p(1).lut = a(1);
end
p(1).lb_step = stepSize(p(1).lb); 
p(1).ub_step = stepSize(p(1).ub);
p(1).alpha = 1; % opaque
p(1).smooth = false;
p(1).interp = 1; % nearest
p(1).volume = 1; % first volume

%% estimate StepSize for java spinner
function d = stepSize(val)
d = abs(val/10);
% d = round(d, 1, 'significant');
d = str2double(sprintf('%.1g', d));
d = max(d, 0.01);
if d>4, d = round(d/2)*2; end

%% Return nii struct from nii struct, nii fname or other convertible files
function nii = get_nii(fname)
if isstruct(fname), nii = fname; return;
elseif iscellstr(fname), nam = fname{1};
else, nam = fname;
end
try nii = nii_tool('load', strtrim(nam));
catch, nii = dicm2nii(fname, pwd, 'no_save');
end

%% Get figure/plot position from FoV for layout
% siz is in pixels, while pos is normalized.
function [siz, pos] = plot_pos(mm, layout)
if layout==1 % 1x3
    siz = [sum(mm([2 1 1]))+mm(1)/4 max(mm(2:3))]; % image area width/height
    y0 = mm(2) / siz(1); % normalized width of sag images
    x0 = mm(1) / siz(1); % normalized width of cor/tra image
    z0 = mm(3) / siz(2); % normalized height of sag/cor images
    y1 = mm(2) / siz(2); % normalized height of tra image
    pos = [0 0 y0 z0;  y0 0 x0 z0;  y0+x0 0 x0 y1;  y0+x0*2 0 mm(1)/4/siz(1) min(z0,y1)];
elseif layout<4
    siz = [sum(mm(1:2)) sum(mm(2:3))]; % image area width/height
    x0 = mm(1) / siz(1); % normalized width of cor/tra images
    y0 = mm(2) / siz(2); % normalized height of tra image
    z0 = mm(3) / siz(2); % normalized height of sag/cor images
    if layout == 2 % 2x2 sag at (1,2)
        pos = [x0 y0 1-x0 z0;  0 y0 x0 z0;  0 0 x0 y0;  x0 0 1-x0 y0];
    else % ==3:      2x2 sag at (1,2)
        pos = [0 y0 1-x0 z0;  1-x0 y0 x0 z0;  1-x0 0 x0 y0;  0 0 1-x0 y0];
    end
else
    error('Unknown layout parameter');
end
siz = siz / max(siz) * 800;
%%