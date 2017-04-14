function ea_chirp(options)


if options.prefs.machine.chirp
    load chirp
    sound(y(1:1000),Fs);
end