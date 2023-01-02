function D=ea_minmax(D)

D(~isnan(D(:)))=D(~isnan(D(:)))-(min(D(~isnan(D(:)))));
D(~isnan(D(:)))=D(~isnan(D(:)))./max(D(~isnan(D(:))));