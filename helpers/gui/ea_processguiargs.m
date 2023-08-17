function ea_processguiargs(handles,varargin_ext)

h=dbstack;
callingfunction=h(4).name;

if ~isempty(varargin_ext)
    if ischar(varargin_ext{1})
        switch varargin_ext{1}
            case 'loadsubs'
                if ~strcmp(callingfunction,'lead_group')
                    ea_load_pts(handles,varargin_ext{2});
                else
                    ea_load_pts_group(handles,varargin_ext{2});
                end
        end
    end
end