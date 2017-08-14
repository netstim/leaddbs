function     ea_setaxisorigin(options,resultfig)


% attempt to set the origin of the axis to the center of the electrode
% reconstruction. This is a highly optional step and only works if a patient is
% selected, thus whole step is in try/end brackets for now.

ax=resultfig.CurrentAxes;
coords=ea_load_reconstruction(options);
mc=mean(cell2mat(cellfun(@mean,coords,'Uniformoutput',0)'),1); % mean coordinate

axis([mc(1)-20 mc(1)+20 mc(2)-20 mc(2)+20 mc(3)-20 mc(3)+20]);
