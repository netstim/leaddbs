function ea_flatview(resultfig)
% small function to generate a more flat-design like view good for printing
% and figures

if ~exist('resultfig','var')
    resultfig=gcf;
end



patches = findall(resultfig,'Type','patch');

[patches.DiffuseStrength]=deal(0.5);
[patches.AmbientStrength]=deal(0.4);

[patches.SpecularStrength]=deal(0.05);
[patches.SpecularExponent]=deal(40);

[patches.SpecularColorReflectance]=deal(0);

