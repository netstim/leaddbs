function parsedInput = ea_processguiargs(handles,varargin_ext)

parsedInput = [];

h=dbstack;
callingfunction=h(4).name;
if ~isempty(varargin_ext)
    if ischar(varargin_ext{1})
        switch varargin_ext{1}
            case 'loadsubs'
                parsedInput = varargin_ext{2};
                if ~strcmp(callingfunction,'lead_group')
                    ea_load_pts(handles, parsedInput);
                else
                    ea_load_pts_group(handles, parsedInput);
                end
        end
    end
end