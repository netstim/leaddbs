function ea_magetcheck_norm_peers(options,peerfolders)
% function that checks maget peers for normalization with an ANTs function
% and does normalize them if not done before.

for peer=1:length(peerfolders)
    poptions=options;
    [poptions.root,poptions.patientname]=fileparts(peerfolders{peer});
    poptions.root=[poptions.root,filesep];
    poptions=ea_assignpretra(poptions);
    % make sure peer has been normalized using ANTs
    if ~ismember(ea_whichnormmethod([peerfolders{peer},filesep]),ea_getantsnormfuns)
        ea_dumpmethod(poptions, 'normm');
        ea_normalize_ants(poptions)
    else
        % make sure peer's anatomy files have been coregistered.
        if ~ea_seemscoregistered(poptions)
            ea_coregpreopmr(poptions);
        end
    end
end
