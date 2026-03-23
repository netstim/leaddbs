function [WC, T] = LeG_autoElecs_chainoptimization(app, Img, MRInfo, ProjSurfRaw, ...
    voxDim, ThrR, ThrHU, TMax, TMin, StartTime)
% Reuses the count-free candidate/fit pipeline but swaps the initial shaft
% selection stage for a beam-search conflict optimizer.

[WC, T] = LeG_autoElecs_countguided(app, Img, MRInfo, ProjSurfRaw, ...
    voxDim, ThrR, ThrHU, TMax, TMin, StartTime, 'chainoptimization');
end
