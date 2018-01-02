function ea_chirp(options)
if ~exist('options','var')
   options.prefs=ea_prefs; 
end
try
    if options.prefs.env.logtime;
        f=fopen([options.root,options.patientname,filesep,'ea_timelog.txt'],'a');
        fprintf(f,[date,': Process took %.2f seconds\n'],toc(options.tic));
        fclose(f);
    end
end
try
    if options.prefs.machine.chirp
        load(fullfile(matlabroot, 'toolbox/matlab/audiovideo/chirp.mat'));
        S = warning('off', 'MATLAB:audiovideo:audioplayer:noAudioOutputDevice');
        sound(y(1:1000),Fs/2);
        warning(S);
    end
end