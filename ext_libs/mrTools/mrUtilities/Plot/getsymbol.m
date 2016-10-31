% getsymbol.m
%
%      usage: getsymbol(symbolnum)
%         by: justin gardner
%       date: 10/18/04
%    purpose: return a different symbol for a different
%             number
%       e.g.: getsymbol(1)
function retval = getsymbol(symbolnum,symbol)

if ((nargin ~= 1) && (nargin ~= 2))
  help getsymbol
  return
end


symbols = 'xov^sd<>ph';
if (nargin == 1)
  retval = symbols(mod(symbolnum,length(symbols))+1);
else
  retval = [symbols(mod(symbolnum,length(symbols))+1) symbol];
end