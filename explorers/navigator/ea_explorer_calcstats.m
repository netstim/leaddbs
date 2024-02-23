function [vals,pvals] = ea_explorer_calcstats(obj,voxelsORfibers,patientselection,Outcome)
switch voxelsORfibers
    case 'voxels'
        myvals = cellfun(@full, obj.results.efield, 'Uni', 0);
    case 'fibers'
        myvals = cellfun(@full, obj.results.(ea_conn2connid(obj.connectome)).(ea_explorer_method2methodid(obj)).fibsval, 'Uni', 0);
end
if ~exist('Ouctome','var')
    Outcome = obj.responsevar;
end
if ~exist('patientselection','var') % patientselection can be supplied directly (in this case, obj.patientselection is ignored), e.g. for cross-validations.
    patientselection=obj.patientselection;
end

if any(isnan(Outcome))
    nanidx=find(isnan(Outcome));
    patientselection(nanidx)=[];
    Outcome(nanidx)=[];
    disp(['Removed ' num2str(numel(nanidx)) ' Participant(s) because their outcome was NaN.'])
end
    
switch obj.statsettings.stimulationmodel
    case 'Electric Field'
        for side = 1:size(myvals,2)
            myvals{1,side}(myvals{1,side}<obj.statsettings.efieldthreshold) = NaN;
        end
    case 'Sigmoid Field'
        for side = 1:size(myvals,2)
            myvals{1,side}(:,:) = ea_SigmoidFromEfield(myvals{1,side}(:,:));
            myvals{1,side}(myvals{1,side}<obj.statsettings.efieldthreshold) = NaN;
        end
    case 'VTA'
        for side = 1:size(myvals,2)
            myvals{1,side}(myvals{1,side}<obj.statsettings.efieldthreshold) = NaN;
            myvals{1,side}(myvals{1,side}>=obj.statsettings.efieldthreshold) = 1;
        end
end

if size(Outcome,2)==1 % 1 entry per patient, not per electrode
    Outcome=[Outcome,Outcome]; % both sides the same;
end
if obj.mirrorsides
    Outcome=[Outcome;Outcome];
end

myvalsgroup=myvals; %refresh myvals

patientselectiongroup=patientselection;

if obj.mirrorsides
    patientselectiongroup=[patientselectiongroup,patientselectiongroup + length(obj.allpatients)];
end


for side=1:numel(myvalsgroup)
    %% remove streamlines/voxels that did not receive enough stimulation
    Nvalsgroup = sum(~isnan(myvalsgroup{side}(:,patientselectiongroup)),2);
    % remove fibers that are not connected to enough VTAs/Efields or connected
    % to too many VTAs (connthreshold slider)
    myvalsgroup{side}(Nvalsgroup<((obj.statsettings.connthreshold/100)*length(patientselectiongroup)),patientselectiongroup)=NaN;
    if strcmp(obj.statsettings.statfamily,'2-Sample Tests') % for tests with two samples we also investigate the size of the unstimulated group
        myvalsgroup{side}(Nvalsgroup>((1-(obj.statsettings.connthreshold/100))*length(patientselectiongroup)),patientselectiongroup)=NaN;
    end

    % init outputvars
    vals{1,side}=nan(size(myvalsgroup{side},1),1);
    pvals{1,side}=nan(size(myvalsgroup{side},1),1);

    nonempty=sum(myvalsgroup{side}(:,patientselectiongroup),2,'omitnan')>0;
    nonemptyidx=find(nonempty);

    valsin=myvalsgroup{side}(nonempty,patientselectiongroup);
    outcomein=Outcome(:,side);

    disp(['Calculating ' obj.statsettings.stattest ' for side ' num2str(side) '...']);

    switch obj.statsettings.stattest
        case 'N-Map'
            [valsout,psout]=ea_explorer_stats_nmap(valsin);
        case 'Mean-Map'            
            [valsout,psout]=ea_explorer_stats_meanmap(valsin,outcomein);
        case '2-Sample T-Test' % two-sample t-tests / OSS-DBS
            [valsout,psout]=ea_explorer_stats_2samplettest(valsin,outcomein);
        case '1-Sample T-Test'            
            [valsout,psout]=ea_explorer_stats_1samplettest(valsin,outcomein,obj.statsettings.H0);
        case 'Wilcoxon Signed-Rank Test'            
            [valsout,psout]=ea_explorer_stats_signedranktest(valsin,outcomein,obj.statsettings.H0);
        case 'Wilcoxon Rank-Sum Test'
            [valsout,psout]=ea_explorer_stats_ranksumtest(valsin,outcomein);
        case '1-Sample Weighted Regression'            
            [valsout,psout]=ea_explorer_stats_1sampleweightedlinreg(valsin,outcomein,obj.statsettings.H0);
        case '2-Sample Weighted Regression'
            [valsout,psout]=ea_explorer_stats_2sampleweightedlinreg(valsin,outcomein);
        case 'Spearman'
            [valsout,psout]=ea_explorer_stats_spearman(valsin,outcomein);
        case 'Pearson'
            [valsout,psout]=ea_explorer_stats_pearson(valsin,outcomein);
        case 'Proportion Test'
            [valsout,psout]=ea_explorer_stats_proportiontest(valsin,outcomein);
    end
    vals{1,side}(nonemptyidx)=valsout;
    pvals{1,side}(nonemptyidx)=psout;
end
