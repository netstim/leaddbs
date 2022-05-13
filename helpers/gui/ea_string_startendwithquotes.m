function has_quotes=ea_string_startendwithquotes(pth)
%check if string has quotes

starts_with_quote=false;
ends_with_quote=false;
if strcmp(pth(1),'"')
    starts_with_quote=true;
    pth(1)=[];
end
if strcmp(pth(end),'"')
    ends_with_quote=true;
    pth(end)=[];
end

has_quotes=starts_with_quote || ends_with_quote;

