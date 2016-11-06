function  ea_processguiargs(handles,varargin_ext)

if ~isempty(varargin_ext)
    if ischar(varargin_ext{1})
    switch varargin_ext{1}
        case 'loadsubs'
            ea_load_pts(handles,varargin_ext{2});
    end
    end
end