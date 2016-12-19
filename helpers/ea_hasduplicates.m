function [has,repeatedValues]=ea_hasduplicates(X)

uniqueX = unique(X);
countOfX = hist(X,uniqueX);
indexToRepeatedValue = (countOfX~=1);
repeatedValues = uniqueX(indexToRepeatedValue);
if ~isempty(repeatedValues)
    has=1;
else
    has=0;
end