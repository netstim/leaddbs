function atlases=ea_genatlastable(options)
% This function reads in atlases in the eAuto/atlases directory and
% generates a table of all available atlas files.
% Atlastypes:   1 ? LH
%               2 ? RH
%               3 ? both hemispheres (2 files present both in lhs and rhs
%               folder
%               4 ? mixed (one file with one cluster on each hemisphere)
%               5 ? midline (one file with one cluster in total)
%
%
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

lhcell=cell(0); rhcell=cell(0); mixedcell=cell(0); midlinecell=cell(0);
delete([options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep,'*_temp.ni*']);
lhatlases=dir([options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep,'*.ni*']);
for i=1:length(lhatlases); 
    lhcell{i}=lhatlases(i).name; 
end
delete([options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep,'*_temp.ni*']);
rhatlases=dir([options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep,'*.ni*']);
for i=1:length(rhatlases); 
    rhcell{i}=rhatlases(i).name; 
end
delete([options.earoot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep,'*_temp.ni*']);
mixedatlases=dir([options.earoot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep,'*.ni*']);
for i=1:length(mixedatlases); 
    mixedcell{i}=mixedatlases(i).name;
end
delete([options.earoot,'atlases',filesep,options.atlasset,filesep,'midline',filesep,'*_temp.ni*']);

midlineatlases=dir([options.earoot,'atlases',filesep,options.atlasset,filesep,'midline',filesep,'*.ni*']);
for i=1:length(midlineatlases); 
    midlinecell{i}=midlineatlases(i).name;
end

% concatenate lh and rh
    todeletelh=[];
    todeleterh=[];
for i=1:length(lhcell)

    [ism, loc]=ismember(lhcell{i},rhcell);
   if ism
       
      todeletelh=[todeletelh,i];
      todeleterh=[todeleterh,loc];
   end
    
end

bothcell=lhcell(todeletelh);
lhcell(todeletelh)=[];
rhcell(todeleterh)=[];

allcell=[lhcell,rhcell,bothcell,mixedcell,midlinecell];
typecell=[repmat(1,1,length(lhcell)),repmat(2,1,length(rhcell)),repmat(3,1,length(bothcell)),repmat(4,1,length(mixedcell)),repmat(5,1,length(midlinecell))];
atlases.names=allcell;
atlases.types=typecell;