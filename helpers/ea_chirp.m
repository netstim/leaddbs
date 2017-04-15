function ea_chirp(options)
if ~exist('options','var')
   options.prefs=ea_prefs; 
end
try
if options.prefs.machine.chirp
    load chirp
    sound(y(1:1000),Fs/2);
end
end