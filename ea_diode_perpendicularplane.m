function Z = ea_diode_perpendicularplane(normvec,p0,X,Y)
    d = -((normvec(1) .* p0(1)) + (normvec(2) .* p0(2)) + (normvec(3) .* p0(3)));
    Z = (-(normvec(1) .* X)-(normvec(2) .* Y) -d) ./ normvec(3);
end




