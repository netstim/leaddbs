function out=ea_arenopoints4side(points)
    %returns true if there are no valid points for this side (if points is empty or if it is all nan)
    %this is helpful for determining if it is a valid side with .coords_mm, .coords_acpc, .trajectory
    %Enrico Opri, 2020
    out=isempty(points) || all(isnan(points(:)));