function varargout=ea_checkmacaque(varargin)

varargout{1}=''; % make default human
varargout{2}=0;
try
    if isfield(varargin{1},'macaquemodus') % options have been entered
        options=varargin{1};

        varargout{1}='';
        if options.macaquemodus
            varargout{1}=['toolbox',filesep,'macaque',filesep];
        else
            varargout{1}='';
        end
    else % handles have been entered
        handles=varargin{1};
        if getappdata(handles.leadfigure,'macaquemodus')
            varargout{1}=['toolbox',filesep,'macaque',filesep];
        else
            varargout{1}='';
        end
    end

    if nargin>1
        if exist([ea_getearoot,'toolbox',filesep,'macaque'],'file')
            varargout{2}=1;
        else
            varargout{2}=0;
            % macaque toolbox not installed. installing:
            success=ea_checkinstall('macaque');
            if success
                varargout{2}=1;
            end
        end
    end
end
