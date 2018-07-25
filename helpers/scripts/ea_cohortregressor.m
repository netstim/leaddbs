function cohreg=ea_cohortregressor(reg)
% converts group indices (n x 1 vector with 1, 2, 3, etc. entries) to a
% n x max(reg) regressor to clean data for cohorts.


cohreg=zeros(size(reg,1),ea_nanmax(reg));

for c=1:size(cohreg,2)
    cohreg(:,c)=reg==c;
end


