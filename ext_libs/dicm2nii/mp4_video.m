classdef mp4_video < VideoWriter
% Create mp4 video from figure window, like VideoWriter.
%   
%  Example for a moving circle:
%    h = rectangle('Position', [0 0 1 1], 'Curvature', 1, 'FaceColor', 'k');
%    axis equal; axis off; xlim([0 10]); ylim([0 10]); % set axis range
%    vw = mp4_video('movingCircle.mp4', 4); % 4 frames per second
%    for x = 0:9 % circle location range
%        h.Position(1:2) = x; % move along diagonal
%        vw.addFrame(); % add current frame to video buffer
%    end
%    vw.save(); % finish the video
%
% No special toolbox is needed for Windows and OSX. For Linux, ffmpeg is needed.
% In case ffmpeg is not available, the uncompressed avi file will be kept.
% 
% This is based on VideoWriter, but restricts the output to the popular mp4
% video only. The major benefits include enabling Linux support for mp4 format,
% supporting Matlab throught ssh/putty connection, and supporting figure with UI
% controls. The rect feature allows all OS to capture part of the figure so to
% exclude UI control if needed.

% 20200323  Wrote it by Xiangrui.Li at gmail.com
  
  properties
    rect; % [left bottom width height] of video in figure. Default full window
  end
  
  methods
    function obj = mp4_video(filename, frameRate, pos)
        % Construct mp4_video object to write video.
        %  vw = mp4_video(filename, frameRate, rect);
        %   filename: result file name (required). .mp4 will be added if needed.
        %   frameRate: frames per second. Default 30.
        %   pos: Position in pixels to extract frame in figure, or a figure
        %   component, for example panel or axes, to capture.
        %       Default [] meeans full figure window.
        profile = 'MPEG-4';
        if ~ismember(profile, {VideoWriter.getProfiles.Name}) % likely Linux
            profile = 'Uncompressed AVI'; % use ffmpeg to compress later
            [err, ~] = system('ffmpeg -version');
            if err
                fprintf(['OS ffmpeg was not installed or not on path. '...
                    'An uncompressed avi video will be created.\n']);
            end
        end
        if ~endsWith(filename, '.mp4'), filename = [filename '.mp4']; end
        obj@VideoWriter(filename, profile); % superclass constructor
        if nargin>2
            if ~isnumeric(pos)
                pos = round(getpixelposition(pos, 1));
                pos(3:4) = floor(pos(3:4)/2) * 2;
            end
            obj.rect = pos;
        end
        if nargin>1 && ~isempty(frameRate), obj.FrameRate = frameRate; end
    end
    
    function addFrame(this, fh)
        % Capture a frame from specified figure or gcf if omitted.
        persistent use_getframe noui;
        if isempty(use_getframe), use_getframe = true; noui = false; end
        if nargin<2 || ~ishandle(fh), fh = gcf; end
        if this.FrameCount < 1 % determine if use getframe and -noui
            open(this);
            try getframe(fh); catch, use_getframe = false; end
            if ~use_getframe
               try print(fh, '-RGBImage', '-r0'); catch, noui = true; end
            end
        end
        
        if use_getframe % print() works for all, but getframe() is faster
            if isempty(this.rect), A = getframe(fh); % faster without rect
            else, A = getframe(fh, this.rect);
            end
        elseif noui % needed when display disabled
            A = print(fh, '-RGBImage', '-r0', '-noui');
        else % try with UI, only use '-noui' if needed
            A = print(fh, '-RGBImage', '-r0');
        end
        this.writeVideo(A);
    end
    
    function save(this)
        % Finish writing and close video file.
        this.close();
        if strcmpi(this.FileFormat, 'mp4'), return; end % suppose getframe works 
        sz = [this.Width this.Height];
        if isempty(this.rect), this.rect = [0 0 sz-mod(sz,2)]; end
        opts = '-y -pix_fmt yuv420p '; % overwrite, yuv420p for "dumb players"
        if ~isequal(this.rect(3:4), sz)
            r4 = this.rect([3 4 1 2]); r4(4) = sz(2)-r4(4)-r4(2);
            opts = sprintf('%s-vf crop=%d:%d:%d:%d ', opts, r4);
        end
        nam = fullfile(this.Path, this.Filename);
        [err, str] = system(['ffmpeg -i "' nam '" ' opts '"' nam(1:end-4) '"']);
        if err, disp(str); else, delete(nam); end
    end
    
    function set.rect(this, rect)
        if isempty(rect), return; end
        rect = floor(rect);
        if any(mod(rect(3:4),2))
            error('Width and height must be even numbers');
        end
        this.rect = rect;
    end
  end
  
  methods(Hidden)
    % Override inherited methods, make them hidden
    function lh = addlistener(varargin); lh=addlistener@VideoWriter(varargin{:}); end
    function lh = listener(varargin); lh=listener@VideoWriter(varargin{:}); end
    function p = findprop(varargin); p = findprop@VideoWriter(varargin{:}); end
    function lh = findobj(varargin); lh = findobj@VideoWriter(varargin{:}); end
    function TF = eq(varargin); TF = eq@VideoWriter(varargin{:}); end
    function TF = ne(varargin); TF = ne@VideoWriter(varargin{:}); end
    function TF = lt(varargin); TF = lt@VideoWriter(varargin{:}); end
    function TF = le(varargin); TF = le@VideoWriter(varargin{:}); end
    function TF = gt(varargin); TF = gt@VideoWriter(varargin{:}); end
    function TF = ge(varargin); TF = ge@VideoWriter(varargin{:}); end
    function lh = get(varargin); lh = get@VideoWriter(varargin{:}); end
    function notify(varargin); notify@VideoWriter(varargin{:}); end
    function delete(varargin); delete@VideoWriter(varargin{:}); end
    function close(varargin); close@VideoWriter(varargin{:}); end
    function open(varargin); open@VideoWriter(varargin{:}); end
    function set(varargin); set@VideoWriter(varargin{:}); end
    function saveobj(varargin); saveobj@VideoWriter(varargin{:}); end
    function setdisp(varargin); setdisp@VideoWriter(varargin{:}); end
    function writeVideo(varargin); writeVideo@VideoWriter(varargin{:}); end
  end
end
