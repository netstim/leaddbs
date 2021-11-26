function D=ea_minmax(D)

D(:)=D(:)-(min(D(:)));
D(:)=D(:)./max(D(:));