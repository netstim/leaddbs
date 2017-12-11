function ea_warning(msg)
% the warning dialog of this function may be silenced by adding
% evalin('base','WARNINGSILENT=1;');
% to the beginning of your code.
warning(msg);
try
    ws=evalin('base','WARNINGSILENT');
    if ws
        return
    end
end
errordlg(sprintf([msg, '\n\n[This warning dialog may be silenced by adding \n\n''evalin(''base'',''WARNINGSILENT=1;'');''\n\n to the beginning of your code.]']));
