function mytests = ea_explorer_stats_list
    warning('off','all')
    mytests = table('Size',[0,5],...
        'VariableTypes',{'categorical','string','categorical','cell','cell'},...
        'VariableNames',{'name','file','type','outcometype','compatibility'});
    row=0;
    %%
    row=row+1;
    mytests.name(row) = "N-Map";
    mytests.file(row) = "ea_explorer_stats_nmap.m";
    mytests.type(row) = "Descriptive";
    mytests.outcometype{row} = {'gradual','binary'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field','VTA'};
    %%
    row=row+1;
    mytests.name(row) = "Mean-Map";
    mytests.file(row) = "ea_explorer_stats_meanmap.m";
    mytests.type(row) = "Descriptive";
    mytests.outcometype{row} = {'gradual','binary'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field','VTA'};    
    %%
    row=row+1;
    mytests.name(row) = "1-Sample T-Test";
    mytests.file(row) = "ea_explorer_stats_1samplettest.m";
    mytests.type(row) = "1-Sample Tests";
    mytests.outcometype{row} = {'gradual','binary'};
    mytests.compatibility{row} = {'VTA'};
    %%
    row=row+1;
    mytests.name(row) = "2-Sample T-Test";
    mytests.file(row) = "ea_explorer_stats_2samplettest.m";
    mytests.type(row) = "2-Sample Tests";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'VTA'};
    %%
    row=row+1;
    mytests.name(row) = "Spearman";
    mytests.file(row) = "ea_explorer_stats_spearman.m";
    mytests.type(row) = "Correlations";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field'};    
    %%
    row=row+1;
    mytests.name(row) = "Pearson";
    mytests.file(row) = "ea_explorer_stats_pearson.m";
    mytests.type(row) = "Correlations";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field'};    
    %%
    row=row+1;
    mytests.name(row) = "1-Sample Weighted Regression";
    mytests.file(row) = "ea_explorer_stats_1sampleweightedlinreg.m";
    mytests.type(row) = "1-Sample Tests";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field','VTA'};
    %%
    row=row+1;
    mytests.name(row) = "2-Sample Weighted Regression";
    mytests.file(row) = "ea_explorer_stats_1sampleweightedlinreg.m";
    mytests.type(row) = "2-Sample Tests";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'Electric Field','Sigmoid Field','VTA'};
    %%
    row=row+1;
    mytests.name(row) = "Wilcoxon Signed-Rank Test";
    mytests.file(row) = "ea_explorer_stats_signedranktest.m";
    mytests.type(row) = "1-Sample Tests";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'VTA'};
    %%
    row=row+1;
    mytests.name(row) = "Wilcoxon Rank-Sum Test";
    mytests.file(row) = "ea_explorer_stats_ranrumtest.m";
    mytests.type(row) = "2-Sample Tests";
    mytests.outcometype{row} = {'gradual'};
    mytests.compatibility{row} = {'VTA'};
end