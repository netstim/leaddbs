function is=ea_isbinary(I)

I=I(~isnan(I));
I=I(~isinf(I));

is=isequal(double(I),double(logical(I)));


