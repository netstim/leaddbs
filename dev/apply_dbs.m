function potential = apply_dbs(vol,elec,val,unipolar,constvol)
if constvol
    if unipolar
        dirinodes = get_surf_nodes(vol.hex);
        dirinodes = [dirinodes, elec];
    else
        dirinodes = elec;
    end
    
    rhs = zeros(length(vol.pos),1);
    dirival = zeros(size(vol.pos,1),1);
    dirival(elec) = val;
else
    dirinodes = 1;
    dirival = zeros(size(vol.pos,1),1);
    
    rhs = zeros(size(vol.pos,1),1);
    
    if unipolar
        catnodes = get_surf_nodes(vol.hex);
        rhs(elec) = val;
        rhs(catnodes) = -sum(val)/length(catnodes);
    else
        rhs(elec) = val;
        if abs(sum(val > 1e-10))
            warning('Sum of current is not zero! This leads to an inaccurate simulation!');
        end
    end
end

[stiff rhs] = dbs(vol.stiff,rhs,dirinodes,dirival);

potential = sb_solve(stiff,rhs);    
end
