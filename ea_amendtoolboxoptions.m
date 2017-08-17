function options=ea_amendtoolboxoptions(options)

if ~isfield(options,'lc') % might be predefined from an exported script..
    try
        options.lc = options.prefs.machine.lc;
    catch
        options.lc = [];
    end
end

% append 2D options.
options.d2 = ea_tdhandles2options([], options.d2);
