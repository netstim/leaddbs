function suffx=ea_getarchsuffx
if ispc
    suffx='exe';
else
    suffx=computer('arch');
end