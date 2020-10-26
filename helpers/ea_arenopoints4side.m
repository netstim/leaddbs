function out=ea_arenopoints4side(coords_list, side)
    %returns true if there are no valid points for this side (if points is empty or if it is all nan)
    %this is helpful for determining if it is a valid side with .coords_mm, .coords_acpc, .trajectory
    %Enrico Opri, 2020
    
    %this assumes hardcoded R side as 1, and L as 2
    if length(coords_list)<side
        %the side does not exist
        out=true;
    else
        out=ea_arenopoints(coords_list{side});
    end
    
    function out=ea_arenopoints(points)
        out=isempty(points) || all(isnan(points(:)));
    end
end