function varargout=ea_checkmacaque(varargin)

options=varargin{1};

varargout{1}='';
try
    if options.macaquemodus
        varargout{1}=['toolbox',filesep,'macaque',filesep];
    else
        varargout{1}='';
    end
end

if nargin>1
    if exist([options.earoot,'toolbox',filesep,'macaque'],'file')
        
     varargout{2}=1;
    else
        varargout{2}=0;
    end
end
        
        
        
        
        