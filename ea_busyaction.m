function ea_busyaction(varargin)
% function displays or hides a spinning wheel on top right corner (default)
% of figure.

if isa(varargin{2}, 'matlab.ui.control.Image')
    switch varargin{1}
        case 'on'
            varargin{2}.Visible = 'on';
            varargin{2}.ImageSource = fullfile(ea_getearoot, 'icons', 'busy.gif');
            varargin{2}.Tooltip = 'Busy';
        case 'off'
            varargin{2}.ImageSource = fullfile(ea_getearoot, 'icons', 'idle.png');
            varargin{2}.Tooltip = 'Idle';
    end
    drawnow;
    return;
end

try
    onoff=varargin{1};
    fighandle=varargin{2};

    if nargin>2
        pos=varargin{3};
    else
        sz=get(fighandle,'Position');
        pos=[sz(3)-80-10 sz(4)-10-80 50 50];
    end

    if ischar(pos)
        switch pos
            case 'dbs'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead DBS (busy...)';
                offfigtit='Lead DBS';
            case 'coreg'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-42-2 50 50];
                onfigtit='Check Coregistration (busy...)';
                offfigtit='Check Coregistration';
            case 'atlcontrol'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50+4 sz(4)-50+8 50 50];
                onfigtit='Atlas Control (busy...)';
                offfigtit='Atlas Control';
            case 'anatomy'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead Anatomy (busy...)';
                offfigtit='Lead Anatomy';
            case 'mapper'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead Connectome Mapper (busy...)';
                offfigtit='Lead Connectome Mapper';
            case 'connectome'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead Connectome (busy...)';
                offfigtit='Lead Connectome';
            case 'predict'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead Predict (busy...)';
                offfigtit='Lead Predict';
            case 'or'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Lead OR (busy...)';
                offfigtit='Lead OR';
            case 'acpc'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='ACPC/MNI-space conversions (busy...)';
                offfigtit='ACPC/MNI-space conversions';
            case 'group'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-5 sz(4)-50 50 50];
                onfigtit='Lead Group Analysis (busy...)';
                offfigtit='Lead Group Analysis';
            case 'reco'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Manual reconstruction (busy...)';
                offfigtit='Manual reconstruction';
            case 'stim'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50 sz(4)-39 50 50];
                onfigtit='Stimulation Parameters (busy...)';
                offfigtit='Stimulation Parameters';
            case 'wavelet'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Wavelet-Based Image Fusion (busy...)';
                offfigtit='Wavelet-Based Image Fusion';
            case 'normcheck'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
            case 'normcheckstructures'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Check registration of specific structures (busy...)';
                offfigtit='Check registration of specific structures';
            case 'trajectory'
                sz=get(fighandle,'Position');
                pos=[sz(3)-50-2 sz(4)-50-2 50 50];
                onfigtit='Edit Trajectory (busy...)';
                offfigtit='Edit Trajectory';
        end
    end

    switch onoff
        case 'on'
            if ~exist('onfigtit','var')
                figtit=get(fighandle,'Name');
                set(fighandle,'name',[figtit,' (busy...)']);
            else
                set(fighandle,'name',onfigtit);
            end

            spinner=getappdata(fighandle,'spinner');

            if isempty(spinner)
                try
                    % R2010a and newer
                    iconsClassName = 'com.mathworks.widgets.BusyAffordance$AffordanceSize';
                    iconsSizeEnums = javaMethod('values',iconsClassName);
                    SIZE_32x32 = iconsSizeEnums(2);  % (1) = 16x16,  (2) = 32x32
                    spinner = com.mathworks.widgets.BusyAffordance(SIZE_32x32, 'Busy');  % icon, label
                catch
                    % R2009b and earlier
                    redColor   = java.awt.Color(1,0,0);
                    blackColor = java.awt.Color(0,0,0);
                    spinner = com.mathworks.widgets.BusyAffordance(redColor, blackColor);
                end
                spinner.getComponent.setBackground(java.awt.Color(1, 1, 1));
                spinner.setPaintsWhenStopped(true);  % default = false
                spinner.useWhiteDots(false);         % default = false (true is good for dark backgrounds)
            end
            ea_javacomponent(spinner.getComponent, pos, fighandle);
            spinner.setBusyText('Busy');
            spinner.start;

            setappdata(fighandle,'spinner',spinner);

            % lock mouse pointer in non-dev environment:
            prefs = ea_prefs;
            if ~prefs.env.dev
                set(fighandle, 'pointer', 'watch');
            end

            drawnow;

        case 'off'
            if ~exist('offfigtit','var')
                figtit=get(fighandle,'Name');
                set(fighandle,'name',strrep(figtit,' (busy...)', ''));
            else
                set(fighandle,'name',offfigtit);
            end

            spinner=getappdata(fighandle,'spinner');
            spinner.stop;
            spinner.setBusyText('Idle');
            %spinner.getComponent.setVisible(false);

            % change mousewheel, too:
            set(fighandle, 'pointer', 'arrow');


            disp('** Process done.');

        case 'del'
            if ~exist('offfigtit','var')
                figtit=get(fighandle,'Name');
                set(fighandle,'name',strrep(figtit,' (busy...)', ''));
            else
                set(fighandle,'name',offfigtit);
            end

            spinner=getappdata(fighandle,'spinner');
            spinner.stop;
            spinner.setBusyText('Idle');
            spinner.getComponent.setBackground(java.awt.Color(0,0,0));
            [hjObj, hContainer] = ea_javacomponent(spinner.getComponent, pos, fighandle);
            delete(hContainer);
            spinner.getComponent.setVisible(false)
            setappdata(fighandle,'spinner',[]);
            % change mousewheel, too:
            set(fighandle, 'pointer', 'arrow');
            disp('** Process done.');
    end
end
