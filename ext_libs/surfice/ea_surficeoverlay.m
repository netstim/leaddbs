function ea_surficeoverlay(overlay, threshold, sideCode, showColorbar, colorbarPosition, useSmoothedMesh)
% Wrapper to export mesh+overlay image from Surf-Ice

% Check overlay path
if ischar(overlay)
    overlay = {overlay};
end
overlay = cellfun(@GetFullPath, overlay, 'Uni', 0);

% Use auto threshold by default
if ~exist('threshold','var') || isempty(threshold)
    autoThresh = 1;
else
    autoThresh = 0;
end

% Determine threshold
if autoThresh
    threshold = ea_sfc_getautothresh(overlay);
elseif size(threshold,1) ~= length(overlay) % Use unified threshold for all overlays
	threshold = repmat(threshold, length(overlay), 1);
end

% Determine sides of the exported images
if ~exist('sideCode','var')
    sideCode = {'r', 'l'}; % Left, Right, Superior
elseif isnumeric(sideCode)
    allSideCode = {'r', 'l', 'cb', 's'}; % Left, Right, Cerebellum, Superior
    sideCode = allSideCode(sideCode);
end

% No colorbar by default
if ~exist('showColorbar','var')
    showColorbar = 0;
end

% Top colorbar by default
if ~exist('colorbarPosition','var')
    colorbarPosition = 4;
end

% Use normal mesh instead of smoothed mesh by default
if ~exist('useSmoothedMesh','var')
    useSmoothedMesh = 0;
end

% Get mesh file path
if ischar(useSmoothedMesh) % Custom mesh supplied
    mesh = useSmoothedMesh;
elseif useSmoothedMesh
    mesh = [ea_space,'surf_smoothed.rh.mz3'];
else
    mesh = [ea_space,'surf.rh.mz3'];
end

% Construct script for Surf-Ice
script = ['import gl;', ...
    'gl.resetdefaults();', ...
    'gl.orientcubevisible(0);', ...
    'gl.colorbarvisible(', num2str(showColorbar), ');',...
    'gl.colorbarposition(', num2str(colorbarPosition), ')', ...
    'gl.meshload(''', mesh, ''');'];

for f=1:length(overlay)
    % First two values in threshold are in positive range
    script = [script, 'gl.overlayload(''', overlay{f}, ''');', ...
        'gl.overlayminmax(1, ', num2str(threshold(f,1)), ', ', num2str(threshold(f,2)), ');', ...
        'gl.overlaycolorname(1, ''Red-Yellow'');'];

    % Last two values in threshold are in negative range
    if numel(threshold(f,:))>2 && ~isnan(threshold(f,3))
        script = [script, 'gl.overlayload(''', overlay{f}, ''');', ...
            'gl.overlayminmax(1, ', num2str(threshold(f,3)), ', ', num2str(threshold(f,4)), ');', ...
            'gl.overlaycolorname(1, ''Blue-Green'');'];
    end

    % Base path of images to be exported (.nii.gz or .nii stripped)
    overlayBasePath = ea_niifileparts(overlay{f});

    % Export images for right side
    if ismember('r', sideCode)
        script = [script, 'gl.meshhemisphere(1);', ...
            'gl.azimuthelevation(90,0);', ...
            'gl.savebmp(''', overlayBasePath, '_r_lat.png'');', ...
            'gl.azimuthelevation(-90,0);', ...
            'gl.savebmp(''', overlayBasePath, '_r_med.png'');'];
    end

    % Export images for left side
    if ismember('l', sideCode)
        script = [script, 'gl.meshhemisphere(-1);', ...
            'gl.azimuthelevation(-90,0);', ...
            'gl.savebmp(''', overlayBasePath, '_l_lat.png'');', ...
            'gl.azimuthelevation(90,0);', ...
            'gl.savebmp(''', overlayBasePath, '_l_med.png'');'];
    end

    % Export images for cerebellum
    if ismember('cb', sideCode)
        script = [script, 'gl.meshhemisphere(0);', ...
            'gl.azimuthelevation(0,-45);', ...
            'gl.savebmp(''', overlayBasePath, '_cb.png'');'];
    end

    % Export images for superior side
    if ismember('s', sideCode)
        script = [script, 'gl.meshhemisphere(0);', ...
            'gl.azimuthelevation(0,90);', ...
            'gl.savebmp(''', overlayBasePath, '_s.png'');'];
    end

    % Close current overly
    script = [script, 'gl.overlaycloseall();'];
end

% Quit after finished
script = [script, 'gl.quit();'];

% Run scritpt using Surf-Ice
ea_surfice(script, 1);

% Wait for image writing
pause(0.5);

% Crop files
for f=1:length(overlay)
    overlayBasePath = ea_niifileparts(overlay{f});
    [overlayPath, overlayName] = fileparts(overlayBasePath);

    overlayImages = ea_regexpdir(overlayPath, [overlayName, '_(r_lat|r_med|l_lat|l_med|cb|s)\.png'], 0);
    for i=1:length(overlayImages)
        [image, transparency] = ea_crop_img(overlayImages{i}, 5);
        imwrite(image, overlayImages{i}, 'Alpha', transparency);
    end
end
