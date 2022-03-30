function ea_warning(msg)
% The warning dialog of this function may be silenced by adding
% evalin('base','WARNINGSILENT=1;'); to the beginning of your code.

msg = sprintf(msg);

warning(msg);

if ismember('WARNINGSILENT', evalin('base', 'who')) && evalin('base','WARNINGSILENT')
    return;
else
    errordlg(msg, '');
end
