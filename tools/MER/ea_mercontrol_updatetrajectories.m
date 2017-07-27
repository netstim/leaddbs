function ea_mercontrol_updatetrajectories(handles, side_str)
if ~exist('side_str', 'var')
    side_str = 'both';
end
[~, side_ids, ~] = ea_detsidestr(side_str);

resultfig = getappdata(handles.mercontrolfig, 'resultfig');
options = getappdata(handles.mercontrolfig, 'options');
merstruct = getappdata(resultfig, 'merstruct');
elspec = getappdata(resultfig, 'elspec');

if options.native
    spc = 'native';
else
    spc = 'mni';
end


for sid = side_ids
    im_depth = merstruct.implant_depth(sid);
    im_mm = merstruct.dbs_coords_mm.(spc){sid};
    
    elec_diff = mean(diff(im_mm));  % Average 3-D distance through electrode contacts.
    elec_uv = elec_diff / norm(elec_diff); % unit vector through electrode contacts if tip was at origin.
    im_mm = bsxfun(@minus, im_mm, elec_uv * elspec.contact_length / 2); %shift coordinates to top of contact.
    
    side_transl = merstruct.translations.(spc){sid};
    
    for pos_ix = 1:length(merstruct.tract_info)
        pos = merstruct.tract_info(pos_ix).label;
        curr_dist = merstruct.currentmer.(pos).dist(sid) - im_depth;
        startpoint = im_mm(1, :) + side_transl(pos_ix, :) + (elec_uv .* curr_dist);
        endpoint = startpoint + elec_uv * merstruct.length;
        stepsize = (endpoint - startpoint) / (merstruct.n_pnts - 1);
        pos = merstruct.tract_info(pos_ix).label;
        merstruct.currentmer.(pos).trajectory{sid} = bsxfun(@plus,...
            (0:merstruct.n_pnts-1)' * stepsize, startpoint);
    end
end
  
setappdata(resultfig,'merstruct',merstruct)