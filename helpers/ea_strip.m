function strOut = ea_strip(strIn)
% EA_STRIP Trims whitespace and specific non-printable characters from start and end.
%
% This function removes:
%   - Standard whitespace and control characters: ASCII 0-32 (includes Space, Tab, LF, CR, NAK)
%   - ASCII 127: Delete (DEL)
%   - ASCII 133: Next Line (NEL)
%   - ASCII 160: Non-breaking space (NBSP)
%   - ASCII 8199: Figure Space
%   - ASCII 8239: Narrow No-Break Space

if isempty(strIn)
    strOut = strIn;
    return;
end

% \x00-\x20       : Range of ASCII control characters and space (0-32)
% \x7F            : ASCII 127 (DEL)
% \x85            : ASCII 133 (NEL)
% \xA0            : ASCII 160 (NBSP)
% \x{2007}        : Unicode 8199 (Figure Space)
% \x{202F}        : Unicode 8239 (Narrow No-Break Space)    
pattern = '^[\x00-\x20\x7F\x85\xA0\x{2007}\x{202F}]+|[\x00-\x20\x7F\x85\xA0\x{2007}\x{202F}]+$';

% regexprep handles char arrays, string arrays, and cell arrays automatically
strOut = regexprep(strIn, pattern, '');
