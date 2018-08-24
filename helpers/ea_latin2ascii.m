function output = ea_latin2ascii(input)
% Convert latin string [cell] to ascii string [cell]
%
% TODO: extend to other special characters

latin = {'ä', 'ö', 'ü', 'Ä', 'Ö', 'Ü', 'ß'};
ascii = {'ae', 'oe', 'ue', 'Ae', 'Oe', 'Ue', 'ss'};

output = regexprep(input, latin, ascii);
