function [numStr, mant, expo]= norm_numbers(num, expOffset, mantLength, sepStr, expoStep)
%function [numStr, mant, expo]= norm_nummbers(num, expOffset, mantLength)

if ~exist('mantLength')
    mantLength= 3;   
end
if ~exist('expOffset')
    expOffset= 1;   
end
if ~exist('expoStep')
    expoStep= 1;
end
if ~exist('sepStr')
    sepStr= '*10^';
end
if num == 0;
    mant= 0;    expo= 0;
else
    expo= floor(floor(log10(abs(num)) - expOffset + 1)/expoStep)*expoStep;
    mant= num*(10^-expo);
end

printStr= strcat('%', num2str(expOffset), '.', num2str(mantLength), 'f', sepStr, '%d');
numStr= sprintf(printStr, mant, expo);