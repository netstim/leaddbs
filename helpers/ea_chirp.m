function ea_chirp(options)
if ~exist('options','var')
   options.prefs=ea_prefs; 
end
try
    if options.prefs.machine.chirp
        load(fullfile(matlabroot, 'toolbox/matlab/audiovideo/chirp.mat'));
        sound(y(1:1000),Fs/2);
    end
end