function ea_unifiedmapping_crossval_visualize(explorer,I,Ihat,cvs,posthoccorrectforgroup,sel,group)

if ~exist('group','var')
    group=nan;
end

% check if groups exist:
if length(unique(explorer.M.patient.group))>1 && posthoccorrectforgroup
    Ihat=ea_resid(ea_cohortregressor(explorer.M.patient.group(sel)),Ihat);
    suffix=' [Corrected for Cohort]';
else
    suffix='';
end

groupID = explorer.M.patient.group(sel);
groupColors = explorer.M.groups.color(unique(groupID),:);
groupColors = rgb2hsv(groupColors);
groupColors(:,2) = groupColors(:,2)+((1-groupColors(:,2))/2);
groupColors = hsv2rgb(groupColors);
empiricallabel = [strrep(explorer.responsevarlabel, '_', ' '), ' (Empirical)'];
if explorer.doactualprediction
    fibscorelabel = [strrep(explorer.responsevarlabel, '_', ' '), ' (Predicted)'];
else
    fibscorelabel=[];
    switch explorer.statmetric
        case {'Correlations / E-fields (Irmen 2020)','Reverse T-Tests / E-Fields (binary vars)'}
            fibscorelabel=[fibscorelabel,'Weighted '];
    end
    switch lower(explorer.basepredictionon)
        case 'mean of scores'
            fibscorelabel=[fibscorelabel,'Mean of'];
        case 'sum of scores'
            fibscorelabel=[fibscorelabel,'Sum of'];
        case 'peak of scores'
            fibscorelabel=[fibscorelabel,'Peak of'];
        case 'peak 5% of scores'
            fibscorelabel=[fibscorelabel,'Peak 5% of'];
        case 'profile of scores: spearman'
            fibscorelabel=[fibscorelabel,'Spatial Spearman Correlation of'];
        case 'profile of scores: pearson'
            fibscorelabel=[fibscorelabel,'Spatial Pearson Correlation of'];
        case 'profile of scores: bend'
            fibscorelabel=[fibscorelabel,'Spatial Bend Correlation of'];
    end
    switch explorer.statmetric
        case {'Two-Sample T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)','One-Sample Tests / VTAs / PAM (OSS-DBS)','Reverse T-Tests / E-Fields (binary vars)'}
            fibscorelabel=[fibscorelabel,' Fiber-T-Scores'];
        case 'Correlations / E-fields (Irmen 2020)'
            fibscorelabel=[fibscorelabel,' Fiber-R-Scores'];
        case 'Proportion Test (Chi-Square) / VTAs (binary vars)'
            fibscorelabel=[fibscorelabel,' Fiber-Chi-Scores'];
        case 'Binomial Tests / VTAs (binary vars)'
            fibscorelabel=[fibscorelabel,' Fiber-Binom-PDF'];
        case 'Plain Connections'
            fibscorelabel=[fibscorelabel,' Plain Connections'];
    end
    fibscorelabel=[fibscorelabel,suffix];
end
if iscell(I)

    for sc=1:length(I)
        fibscorelabel = [strrep(explorer.subscore.labels{sc}, '_', ' '), ' (Predicted)'];
        empiricallabel=[strrep(explorer.subscore.labels{sc}, '_', ' '), ' (Empirical)'];
        % h=ea_corrbox(I{sc},Ihat{sc},0,{['Disc. Fiber prediction ',upper(cvs)],empiricallabel,fibscorelabel},groupID,[],groupColors);
        

        h=ea_corrbox(I{sc},Ihat{sc},0,{['Disc. Fiber prediction ',upper(cvs)],empiricallabel,fibscorelabel});
        try saveas(h,[fileparts(explorer.leadgroup),filesep,'fiberfiltering',filesep,explorer.ID,'_',explorer.subscore.labels{sc},'_',cvs,'_',num2str(explorer.numpcs),'_factor_pca.png']); end
    end

else
    if isnan(group)
        title=['Disc. Fiber prediction ',upper(cvs)];
        groupsuffx='';
    else
        title=['Group ',num2str(group),': Disc. Fiber prediction ',upper(cvs)];
        groupsuffx=['_group_',num2str(group)];
    end
    if exist('pperm', 'var')
        h=ea_corrbox(I(~isnan(Ihat)),Ihat(~isnan(Ihat)),pperm,{title,empiricallabel,fibscorelabel},groupID(~isnan(Ihat)),[],groupColors);
    else
        h=ea_corrbox(I(~isnan(Ihat)),Ihat(~isnan(Ihat)),'permutation',{title,empiricallabel,fibscorelabel},groupID(~isnan(Ihat)),[],groupColors);
    end
    assignin('base','Empirical',I);
    assignin('base','Estimate',Ihat);
    try saveas(h,[fileparts(explorer.leadgroup),filesep,'fiberfiltering',filesep,explorer.ID,'_',explorer.responsevarlabel,'_',cvs,groupsuffx,'.png']); end
end