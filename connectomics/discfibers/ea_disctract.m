classdef ea_disctract < handle
    % Discriminative fiber class to handle visualizations of discriminative fibers in lead dbs resultfig / 3D Matlab figures
    % A. Horn
    
    properties (SetObservable)
        app % app handle to appdesigner view (needs to be supplied to construct)
        ispositive % boolean to distinguish positive or negative fibers (needs to be supplied to construct)
        M % content of lead group project (needs to be supplied to construct)
        tractset % parent object of ea_disctractset
        resultfig % figure handle to plot results
        ID % name / ID of discriminative fibers object
        visible % object visible
        showamount % threshold up to which fibers are shown
        showamountleft
        showamountright
        connthreshold
        statmetric % entry from discfiber settings as initially specified in prefs.machine.lg (which are separately stored for each analysis/object).
        cmap % colormap for fiber
        color % main color
        % the  main content variables:
        
        results 
        % includes subfields by results.connectomename.ttests /
        % results.connectomename.efields with
        % fibcell % cell of all fibers connected
        % fibsval % connection weight value for each fiber to each VTA
        % fibweights % usually T- or R-values associated with each tract
        %
        patientselection % selected patients to include. Note that connected fibers are always sampled from all (& mirrored) VTAs of the lead group file
        allpatients % list of all patients (as from M.patient.list)
        mirrorsides % flag to mirror VTAs / Efields to contralateral sides using ea_flip_lr_nonlinear()
        variables % redundancy protocol only, first entry is target variable, further are covariates
        varlabels % redundancy protocol only, first entry is target variable, further are covariates
        leadgroup % redundancy protocol only, path to original lead group project
        connectome % redundancy protocol only, name of underlying connectome
    end
    
    properties (Access = private)
        switchedFromSpace=3 % if switching space, this will protocol where from
    end
    
    methods
        function obj=ea_disctract(pobj) % class constructor, pobj needs to have M, ispositive and app entries. If it has more, it is treated as a loaded object with all entries.
            if isfield(pobj,'fibcell')
                obj=pobj; % e.g. loaded object from disk supplied.
            else
                % set up data from pobj.M and pobj.app
                obj.M=pobj.M;
                obj.app=pobj.app;
                obj.ispositive=pobj.ispositive;
                obj.resultfig = gcf; % change
                obj.ID = 'mytract'; % change
                if obj.ispositive
                    obj.visible = obj.app.ShowPositiveFibersCheckBox.Value;
                    obj.showamount = obj.app.ShowAmountSliderPositive.Value;
                    obj.showamountleft = obj.app.LeftSliderPositive.Value;
                    obj.showamountright = obj.app.RightSliderPositive.Value;
                    obj.color = [1,0,0];
                    obj.cmap = ea_colorgradient(1024, [1,1,1], obj.color); % white to red, default for positives
                else
                    obj.visible = obj.app.ShowNegativeFibersCheckBox.Value;
                    obj.showamount = obj.app.ShowAmountSliderNegative.Value;
                    obj.showamountleft = obj.app.LeftSliderNegative.Value;
                    obj.showamountright = obj.app.RightSliderNegative.Value;
                    obj.color = [0,0,1];
                    obj.cmap = ea_colorgradient(1024, obj.color, [1,1,1]); % blue to white, default for negatives
                end
                obj.patientselection = obj.M.ui.listselect;
                obj.allpatients = obj.M.patient.list;
                obj.varlabels = obj.M.clinical.labels;
                obj.connectome = obj.M.ui.connectomename;
                obj.connthreshold = obj.app.ConnthresholdSlider.Value;
                
                % set up variables:
                [~,ix]=ismember(obj.app.VariableofInterestDropDown.Value,obj.M.clinical.labels);
                obj.variables = [obj.M.clinical.vars{ix}];

                [hascovariates,covix]=ismember(obj.app.CovariatesListBox.Value,obj.M.clinical.labels);
                if hascovariates
                    obj.variables=[obj.variables,obj.M.clinical.vars{covix}];
                end
                % Check whether Amplitude or VTA volumes should be
                % included, as well.
                if ismember('Stimulation Amplitude',obj.app.CovariatesListBox.Value)
                   warning('This needs to be included.'); 
                end
                
                if ismember('VTA Volume',obj.app.CovariatesListBox.Value)
                    warning('This needs to be included.');
                end
                
                switch obj.app.ModelSetupDropDown.Value
                    case 'T-Tests / VTAs'
                        obj.statmetric=1;
                    case 'Spearman''s Correlation / Efields'
                        obj.statmetric=2;
                end
                
                [fibsweighted,fibsin,obj.fibsval,iaix]=ea_calcdisctract(obj);
                % reformat to fibcell from that:
                
                disp('Reformatting fibers...');
                if(all(diff(fibsweighted(:,4))>=0)) % Fiber indices should be monotonic, reformat using fast method
                    [~,fibiaxfirst]=unique(fibsweighted(:,4),'first');
                    [~,fibiaxlast]=unique(fibsweighted(:,4),'last');
                    fiblen = fibiaxlast - fibiaxfirst + 1;
                    obj.fibcell = mat2cell(fibsweighted(:,1:3),fiblen);
                    obj.fibweights = mat2cell(fibsweighted(:,5),fiblen);
                    obj.fibweights = cellfun(@mean, obj.fibweights);
                else % Fall back to old method if fiber indices are not monotonic
                    fibidx=unique(fibsweighted(:,4));
                    obj.fibcell=cell(length(fibidx),1);
                    obj.fibweights=zeros(length(fibidx),1);
                    cnt=1;
                    ea_dispercent(0,'reformatting fibers');
                    for fib=fibidx'
                        obj.fibcell{cnt}=fibsweighted(fibsweighted(:,4)==fib,1:3);
                        obj.fibweights(cnt)=mean(fibsweighted(fibsweighted(:,4)==fib,5));
                        cnt=cnt+1;
                        ea_dispercent(cnt/length(fibidx));
                    end
                    ea_dispercent(1,'end');
                end
                
                
            end
        end
        
        
        
        function draw(M,discfiberssetting,resultfig,fibsweighted)
            
            I=M.clinical.vars{M.ui.clinicallist}(M.ui.listselect);
            
            % Get discriminative fiber setting
            statmetric = discfiberssetting.statmetric;
            
            % protocol selection to be able to check if same analysis has been run
            % before.
            %opts.percent=connthreshold; % need not protocol anymore, is now dynamic
            %option.
            opts.patientselection=M.ui.listselect;
            opts.regressor=I;
            opts.connectome=M.ui.connectomename;
            opts.allpatients=M.patient.list;
            opts.mirrorsides=M.ui.mirrorsides;
            opts.statmetric=statmetric;
            
            patlist=M.patient.list(M.ui.listselect);
            showfibersset = discfiberssetting.showfibersset;
            pospredthreshold = discfiberssetting.pospredthreshold/100;
            negpredthreshold = discfiberssetting.negpredthreshold/100;
            
            [reforce,connectomechanged,reformat]=ea_discfibers_checkpresence(M,opts); % only static opts need to be equal.
            switch statmetric
                case 1 % ttests
                    savesuffix='_ttests';
                case 2 % spearmans R
                    savesuffix='_spearmansrho';
            end
            if M.ui.mirrorsides
                msuffix='_mirrored';
            else
                msuffix='';
            end
            
            
            
            set(0,'CurrentFigure',resultfig);
            
            % Normalize vals
            vals(isnan(vals))=0;
            % vals=vals./max(abs(vals));
            
            % vals and fibcell to be trimmed for visualization
            tvals=vals;
            tfibcell=fibcell;
            
            % Calculate positive/negative threshold for positive/negative predictive
            % fibers according to 'predthreshold'
            posits=tvals(tvals>0);
            negits=tvals(tvals<0);
            posits=sort(posits,'descend');
            negits=sort(negits,'ascend');
            posthresh=posits(round(length(posits)*pospredthreshold));
            negthresh=negits(round(length(negits)*negpredthreshold));
            
            discfiberID = '';
            
            % Save the original values for reusing in slider
            setappdata(resultfig, ['vals',discfiberID], vals);
            setappdata(resultfig, ['fibcell',discfiberID], fibcell);
            setappdata(resultfig, ['showfibersset',discfiberID], showfibersset);
            setappdata(resultfig, ['pospredthreshold',discfiberID], pospredthreshold);
            setappdata(resultfig, ['negpredthreshold',discfiberID], negpredthreshold);
            setappdata(resultfig, ['posits',discfiberID], posits);
            setappdata(resultfig, ['negits',discfiberID], negits);
            
            switch showfibersset
                case 'positive'
                    negthresh = negits(1)-eps;
                    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ')']);
                case 'negative'
                    posthresh = posits(1)+eps;
                    disp(['Fiber colors: Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh), ')']);
                case 'both'
                    disp(['Fiber colors: Positive (T = ',num2str(posthresh),' ~ ',num2str(posits(1)), ...
                        '); Negative (T = ',num2str(negits(1)),' ~ ',num2str(negthresh),').']);
            end
            
            % Remove tvals and fibers outside the thresholding range
            remove=logical(logical(tvals<posthresh) .* logical(tvals>negthresh));
            tvals(remove)=[];
            tfibcell(remove)=[];
            
            % Rescale positive/negative tvals to [0 1]/[-1 0]
            tvalsRescale = tvals;
            tvalsRescale(tvals>0)=ea_rescale(tvals(tvals>0), [0 1]);
            tvalsRescale(tvals<0)=ea_rescale(tvals(tvals<0), [-1 0]);
            
            % Contruct default blue to red colormap
            colormap(gray);
            fibcmap = ea_colorgradient(1024, [0,0,1], [1,1,1], [1,0,0]);
            setappdata(resultfig, ['fibcmap',discfiberID], fibcmap);
            
            fibcolorInd=tvalsRescale*(size(fibcmap,1)/2-0.5);
            fibcolorInd=fibcolorInd+(size(fibcmap,1)/2+0.5);
            
            % Set alphas of fibers with light color to 0
            colorbarThreshold = 0.60; % Percentage of the pos/neg color to be kept
            negUpperBound=ceil(size(fibcmap,1)/2*colorbarThreshold);
            poslowerBound=floor((size(fibcmap,1)-size(fibcmap,1)/2*colorbarThreshold));
            alphas=zeros(size(fibcolorInd,1),1);
            switch showfibersset
                case 'positive'
                    alphas(round(fibcolorInd)>=poslowerBound) = 1;
                case 'negative'
                    alphas(round(fibcolorInd)<=negUpperBound) = 1;
                case 'both'
                    alphas(round(fibcolorInd)>=poslowerBound) = 1;
                    alphas(round(fibcolorInd)<=negUpperBound) = 1;
            end
            
            alphas(round(fibcolorInd)>=poslowerBound) = 1;
            
            % only show right side
            % for i=1:length(tfibcell)
            %     if any(tfibcell{i}(:,1)<0)
            %         alphas(i) = 0;
            %     end
            % end
            
            fibalpha=mat2cell(alphas,ones(size(fibcolorInd,1),1));
            
            % Plot fibers
            h=streamtube(tfibcell,0.2);
            nones=repmat({'none'},size(fibcolorInd));
            [h.EdgeColor]=nones{:};
            
            % Calulate fiber colors
            colors=fibcmap(round(fibcolorInd),:);
            fibcolor=mat2cell(colors,ones(size(fibcolorInd)));
            
            % Set fiber colors and alphas
            [h.FaceColor]=fibcolor{:};
            [h.FaceAlpha]=fibalpha{:};
            
            % Set colorbar tick positions and labels
            cbvals = tvals(logical(alphas));
            % cbvals=tvalsRescale(logical(alphas));
            switch showfibersset
                case 'positive'
                    cbmap = fibcmap(ceil(length(fibcmap)/2+0.5):end,:);
                    tick = [poslowerBound, length(fibcmap)] - floor(length(fibcmap)/2) ;
                    poscbvals = sort(cbvals(cbvals>0));
                    ticklabel = [poscbvals(1), poscbvals(end)];
                    ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
                case 'negative'
                    cbmap = fibcmap(1:floor(length(fibcmap)/2-0.5),:);
                    tick = [1, negUpperBound];
                    negcbvals = sort(cbvals(cbvals<0));
                    ticklabel = [negcbvals(1), negcbvals(end)];
                    ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
                case 'both'
                    cbmap = fibcmap;
                    tick = [1, negUpperBound, poslowerBound, length(fibcmap)];
                    poscbvals = sort(cbvals(cbvals>0));
                    negcbvals = sort(cbvals(cbvals<0));
                    ticklabel = [min(cbvals), negcbvals(end), poscbvals(1), max(cbvals)];
                    ticklabel = arrayfun(@(x) num2str(x,'%.2f'), ticklabel, 'Uni', 0);
            end
            
            figTitle = 'discfibers';
            discfibersname = ['discfibers', discfiberID];
            cbfigname = ['cbfig', discfiberID];
            discfiberscontrolname = ['discfiberscontrol', discfiberID];
            
            % Plot colorbar
            cbfig = ea_plot_colorbar(cbmap, [], 'h', '', tick, ticklabel);
            set(cbfig, 'NumberTitle', 'off',  'Name', ['Colorbar: ', figTitle]);
            
            % Discriminative fiber control
            discfiberscontrol = ea_discfibers_control(resultfig, discfiberID);
            set(discfiberscontrol, 'NumberTitle', 'off', 'Name', ['Control: ', figTitle]);
            setappdata(discfiberscontrol, 'discfiberID', discfiberID);
            
            setappdata(resultfig, discfibersname, h);
            setappdata(resultfig, cbfigname, cbfig);
            setappdata(resultfig, discfiberscontrolname, discfiberscontrol);
            set(0,'CurrentFigure',resultfig)
        end
        
        
        
    end
    
    methods (Static)
        function changeevent(~,event)
            update_trajectory(event.AffectedObject,event.Source.Name);
        end
    end
end


function [fibsweighted,fibsin,fibsval,iaix]=ea_calcdisctract(obj)

% check that this has not been calculated before:
if ~isempty(obj.results) % something has been calculated
    if isfield(obj.results,conn2connid(obj.connectome))
        
        if isfield(obj.results.(conn2connid(obj.connectome)),method2methodid(obj.statmetric)) % this combination was already calculated.
            return
        end
    end
end

efieldthresh=obj.connthreshold;


% Get discriminative fiber setting
connthreshold = obj.connthreshold/100;
statmetric = obj.statmetric;

options.native = 0;

    allroilist=cell(length(obj.allpatients)*2,2);
    switch statmetric
        case 1 % use paired T-Tests and binary VTA
            suffix='';
        case 2 % use Spearman Rs and E-Fields
            suffix='_efield';
            prefs=ea_prefs;
            if strcmp(prefs.lcm.vatseed,'efield_gauss')
                suffix='_efield_gauss';
            end
    end
    cnt=1;
    for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
        allroilist(cnt,:)={[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_right.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_left.nii']};
        cnt=cnt+1;
    end
    for sub=1:length(obj.allpatients) % all patients - for connected fibers selection ? and always flip
        ea_genflippedjointnii([obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_right.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'vat',suffix,'_left.nii']);
        allroilist(cnt,:)={[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat',suffix,'_left.nii'],[obj.allpatients{sub},filesep,'stimulations',filesep,ea_nt(options),'gs_',obj.M.guid,filesep,'fl_','vat',suffix,'_right.nii']};
        cnt=cnt+1;
    end
    cfile=[ea_getconnectomebase('dMRI'),obj.connectome,filesep,'data.mat'];
    mirroredpatselection=[obj.patientselection,obj.patientselection+length(obj.allpatients)];
    switch statmetric
        case 1 % ttests
            [fibcell,fibsval,XYZmm,nii]=ea_discfibers_heatfibertracts(cfile,{allroilist},mirroredpatselection,{obj.variables(:,1)},connthreshold);
        case 2 % spearmans R
            [fibcell,fibsval,XYZmm,nii,valsmm]=ea_discfibers_heatfibertracts_corr(cfile,{allroilist},mirroredpatselection,{obj.variables(:,1)},efieldthresh);
    end


    
    % Main output of results - this is all we will ever need if statmetric
    % and connectome wont change
    obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibsval=fibsval;
    obj.results.(conn2connid(obj.connectome)).(method2methodid(obj.statmetric)).fibcell=fibcell;

end

function id=method2methodid(method)
switch method
    case 1
        id='ttests';
    case 2
        id='spearman';
end
end

function conname=conn2connid(conname)

conname=strrep(conname,' ','');
conname=strrep(conname,'_','');
conname=strrep(conname,'(','');
conname=strrep(conname,')','');
conname=strrep(conname,'-','');
end



