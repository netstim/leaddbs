function Ihat=ea_ihat_databaselookup_netmap(obj,vals,connval,patientsel,test,training,usemask)

model=connval(patientsel(test),usemask)';
Imodel=obj.responsevar(patientsel(test));

sims=ea_corr(model,vals{1}(:,usemask)',resolvecorrtype(obj.basepredictionon));

sh_exp_sim=7; % exponent by which to sharpen similarities weight
sh_exp_cert=2; % exponent by which to sharpen certainties weight (if available)

% if we don't sharpen (e.g. sh_exp=1) then most patients will be taken into
% account similarly

%% similarities
w1=sims-ea_nanmin(sims(:)); % make all positive
w1=w1.^sh_exp_sim; % sharpen
w1=w1./repmat(sum(w1,2),1,sum(training)); % make sure they sum up to 1

if ~strcmp(obj.certainvarlabel,'None') % certainty sidecar
    %% certainties
    w2=obj.certainvar(patientsel(training));
    w2=repmat(w2,1,sum(test))';
    w2=w2-ea_nanmin(w2(:)); % make all positive
    w2=w2.^sh_exp_cert; % sharpen
    w2=w2./repmat(sum(w2,2),1,sum(training)); % make sure they sum up to 1

    w=w1.*w2;
    w=w./repmat(sum(w,2),1,sum(training));
else

    w=w1;
end

Ihat=w*obj.responsevar(patientsel(training));



function c=resolvecorrtype(base)

switch lower(base)
    case 'spatial correlations (spearman)'
        c='Spearman';
    case 'spatial correlations (pearson)'
        c='Pearson';
    case 'spatial correlations (bend)'
        c='Bend';
end
