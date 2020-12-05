%script: ea_run_extra_machinespecific_init.m
%added by Enrico Opri, 2020.
%
%this script is called at the opening of lead.m method
%here are executed specific initializations that are platform specific
%platform specific init included:
% -prefs.platform.unix.load_shipped_libstdcpp6 (default=false)

%examples of variables available at this point:
% -hObject
% -earoot

%by default this library is unneeded
if prefs.platform.unix.load_shipped_libstdcpp6
    % Include this line if libstdc++.so.6 is needed.
    % However it is preferrable to fix it at system level (e.g. using package build-essentials).
    % Install the matlab-support package and choose to use the system libraries for gcc.
    addpath(fullfile(earoot,'ext_libs/support/unix'));
else
    rmpath(fullfile(earoot,'ext_libs/support/unix'));
end