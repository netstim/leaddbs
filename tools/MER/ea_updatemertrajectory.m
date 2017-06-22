function ea_updatemertrajectory(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
elspec = getappdata(resultfig, 'elspec');

for sid = side_ids
    side_str = side_strs{sid};
    im_depth = merstruct.implant_depth(sid);
    im_mm = merstruct.dbs_contacts_mm{sid};
    
    dxyz = norm(diff(im_mm(1:2, :)));  % total distance between contacts 1 and 2
    slope = mean(diff(im_mm)) / dxyz; %mean(diff(coords_mm{side}))/norm(mean(diff(coords_mm{side})))
    im_mm = bsxfun(@minus, im_mm, slope * elspec.contact_length / 2); %shift coordinates to top of contact.
    
    % TODO: Previous code had check for ea_getnativemni,
    % but result was identical whether == 1 or == 2
    
    side_offs = merstruct.offset * cat(1, merstruct.tract_info.transform);
    if strcmpi(side_str, 'left')
        side_offs(:, 1) = -side_offs(:, 1);
    end
    
    for pos_ix = 1:length(merstruct.tract_info)
        pos = merstruct.tract_info(pos_ix).label;
        offset_implant_mm = bsxfun(@plus, im_mm, side_offs(pos_ix, :));
        curr_dist = merstruct.currentmer.(pos).dist(sid) - im_depth;
        startpoint = offset_implant_mm(1,:) + slope .* curr_dist;
        endpoint = startpoint + slope * merstruct.length;
        stepsize = (endpoint - startpoint) / merstruct.n_pnts;
        merstruct.currentmer.(pos).trajectory{sid} = bsxfun(@plus,...
            (1:merstruct.n_pnts)' * stepsize, startpoint);
    end
    
end
setappdata(resultfig,'merstruct',merstruct)