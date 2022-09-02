function rad=ea_species_adjustsize(rad)
sd=load([ea_space,'spacedef.mat']);

if isfield(sd.spacedef,'species')
    load([ea_getearoot,'common',filesep,'ea_species.mat']); % https://faculty.washington.edu/chudler/facts.html
    humanweight=species{1,2};
    for s=1:size(species,1)
        if iscell(species{s,1})
            scell=species{s,1};
        else
            scell={species{s,1}};
        end
        if ismember(sd.spacedef.species,scell)
            break
        end
    end
    if s==size(species,1) && ~ismember(sd.spacedef.species)
        ea_warning(['Specified species ',sd.spacedef.species,' not found in database.']);
    end
    factor=species{s,2}/humanweight;
else
    factor=1;
end

rad=ceil(rad*(factor));
