function surf_nodes = get_surf_nodes(cell)
connectivity = hist(cell(:),unique(cell));
surf_nodes = find(connectivity ~= 8);
end
