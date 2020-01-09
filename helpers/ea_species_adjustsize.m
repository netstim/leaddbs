function rad=ea_species_adjustsize(rad)
s=load([ea_space,'ea_space_def.mat']);

if isfield(s.spacedef,'species')
    load([ea_getearoot,'common',filesep,'ea_species.mat']); % https://faculty.washington.edu/chudler/facts.html
    humanweight=species{1,2};
    for s=1:size(species,1)
        if iscell(species{s,1})
            scell=species{s,1};
        else
            scell={species{s,1}};
        end
        if ismember(s.spacedef.species,scell)
            break
        end
    end
    if s==size(species,1) && ~ismember(s.spacedef.species)
        ea_warning(['Specified species ',s.spacedef.species,' not found in database.']);
    end
    factor=species{s,2}/humanweight;
else
    factor=1;
end

rad=ceil(rad*(factor));