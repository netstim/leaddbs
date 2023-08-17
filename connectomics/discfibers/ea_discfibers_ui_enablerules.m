function results = ea_discfibers_ui_enablerules(app,flag)
            app.isbusy(1);
            if ~isfield(app.tractset.results,ea_conn2connid(app.tractset.connectome))
                val='off';
            else
                val=ea_bool2onoff(isfield(app.tractset.results.(ea_conn2connid(app.tractset.connectome)),ea_method2methodid(app.tractset)));
            end
            try
                numsides=numel(app.tractset.results.(ea_conn2connid(app.tractset.connectome)).(ea_method2methodid(app.tractset)).fibsval);
            catch % default
                numsides=2;
            end
            % enable controls:
            app.ModelSetupDropDown.Enable = val;
            app.ModelSetupDropDownLabel.Enable = val;


            app.MirrorSidesCheckBox.Enable = val;
            app.VariableofInterestDropDown.Enable = val;

            groupsassigned=(length(unique(app.tractset.M.patient.group))>1);
            if groupsassigned
                app.SplitColorByGroupButton.Enable = val;
            else
                app.SplitColorByGroupButton.Enable = 'off';
                app.SplitColorByGroupButton.Value = 0;
            end

            app.SplitColorBySubscoreButton.Enable = val;
            app.SingleTractAnalysisButton.Enable = val;
            app.SplitColorByPCAButton.Enable = val;


            app.ShowPositiveFibersCheckBox.Enable = val;
            app.ShowNegativeFibersCheckBox.Enable = val;

            app.ColorButtonPositive.Enable = val;
            %%%%%%%%todoapp.ShowNegativeFibersCheckBox.Enable = val;

            app.ConnthresholdSlider.Enable = val;
            app.ConnthresholdLabel.Enable = val;

            app.VariableofInterestDropDownLabel.Enable = val;
            app.CovariatesCheckBox.Enable = val;
            app.ColorBar.Visible = val;
            if strcmp(val, 'off')
                cla(app.ColorBar);
            end
            app.ExportasAtlasButton.Enable = val;
            cvEnable = val;
%             if app.SplitColorByGroupButton.Value
%                 app.ColorBar.Visible = 'off';
%                 cla(app.ColorBar);
%                 app.ExportasAtlasButton.Enable = 'off';
%                 cvEnable = 'off';
%             else
                
            % end

            app.SubscoresListBox.Enable = val;
            app.MultitractWeightsButton.Enable = val;
            app.OpeninCleartuneButton.Enable = val;


            app.ROISwitch.Enable=val;
            app.ROISwitchLabel.Enable=val;
            app.ConnectedFibersSwitch.Enable=val;
            app.ConnectedFibersLabel.Enable=val;

            % CV Panel:
            app.StrategyDropDown.Enable = cvEnable;
            app.StrategyDropDownLabel.Enable = cvEnable;
            app.kEditField.Enable = cvEnable;
            app.kEditFieldLabel.Enable = cvEnable;
            app.RunCVButton.Enable = cvEnable;
            app.PosthoccorrectforgroupCheckBox.Enable = cvEnable;
            app.ModelNormalizationDropDown.Enable = cvEnable;
            app.ModelNormalizationDropDownLabel.Enable = cvEnable;


            switch app.BasePredictiononDropDown.Value
                case {'Histogram','z-scored Histogram'}
                    app.BinsEditField.Enable = cvEnable;
                    app.BinsEditFieldLabel.Enable = cvEnable;
                otherwise
                    app.BinsEditField.Enable = 'off';
                    app.BinsEditFieldLabel.Enable = 'off';
            end

            switch app.StrategyDropDown.Value
                case 'k-fold (randomized)'
                    app.kEditField.Visible = 'on';
                    app.kEditFieldLabel.Visible = 'on';
                    app.kEditFieldLabel.Text = 'k:';
                    app.kEditField.Tooltip = '';
                    app.kEditField.Value = num2str(app.tractset.kfold);
                case 'Leave-Nothing-Out (Permutation-Based)'
                    app.kEditField.Visible = 'on';
                    app.kEditFieldLabel.Visible = 'on';
                    app.kEditFieldLabel.Text = '# perms:';
                    app.kEditField.Tooltip = '';
                    app.kEditField.Value = num2str(app.tractset.Nperm);
                case 'Custom (Sets)'
                    app.kEditField.Visible = 'on';
                    app.kEditFieldLabel.Visible = 'on';
                    app.kEditFieldLabel.Text = '# sets:';
                    app.kEditField.Tooltip = 'Randomly split the data into N sets';
                    app.kEditField.Value = num2str(app.tractset.Nsets);
                otherwise
                    app.kEditField.Visible = 'off';
                    app.kEditFieldLabel.Visible = 'off';
                    app.kEditField.Tooltip = '';
            end

            switch app.StrategyDropDown.Value
                case 'Custom (Patients)'
                    [~, patNames] = fileparts(app.tractset.M.patient.list);
                    app.TrainOnListBox.Items = patNames;
                    app.PredictListBox.Items = patNames;
                case 'Custom (Cohorts)'
                    if length(app.tractset.M.groups.group)>1
                        app.TrainOnListBox.Items=sprintfc('%d',(app.tractset.M.groups.group));
                        app.PredictListBox.Items=sprintfc('%d',(app.tractset.M.groups.group));
                    end
                case 'Import Model'
                    if ~strcmp(app.tractset.ExternalModelFile, 'None')
                        [~, ExternalModelFileName] = fileparts(app.tractset.ExternalModelFile);
                        app.ImportedModelName.Value = ExternalModelFileName;
                    end

                    [~, patNames] = fileparts(app.tractset.M.patient.list);
                    app.PredictListBox.Items = patNames;                    
                case 'Custom (Sets)'
                    if str2double(app.kEditField.Value)>1
                        app.TrainOnListBox.Items=sprintfc('%d',(1:str2double(app.kEditField.Value))');
                        app.PredictListBox.Items=sprintfc('%d',(1:str2double(app.kEditField.Value))');
                    end
                case 'Custom (Subcohorts)'
                    app.TrainOnListBox.Items=app.tractset.setlabels;
                    app.PredictListBox.Items=app.tractset.setlabels;
            end

            switch app.StrategyDropDown.Value
                case 'Custom (Patients)'
                    cvlEnable=cvEnable;
                    if isempty(app.TrainOnListBox.Value) || isempty(app.PredictListBox.Value)
                        app.RunCVButton.Enable = 'off';
                    end
                case 'Import Model'
                    cvlEnable=cvEnable;
                    if isempty(app.PredictListBox.Value)
                        app.RunCVButton.Enable = 'off';
                    end
                case 'Custom (Cohorts)'
                    if length(app.tractset.M.groups.group)>1
                        cvlEnable=cvEnable;
                        if isempty(app.TrainOnListBox.Value) || isempty(app.PredictListBox.Value)
                            app.RunCVButton.Enable = 'off';
                        end
                    else
                        app.RunCVButton.Enable = 'off';
                        cvlEnable='off';
                    end

                case 'Custom (Sets)'
                    if str2double(app.kEditField.Value)>1
                        if isempty(app.TrainOnListBox.Value) || isempty(app.PredictListBox.Value)
                            app.RunCVButton.Enable = 'off';
                        end
                        cvlEnable=cvEnable;
                    else
                        app.RunCVButton.Enable = 'off';
                        cvlEnable='off';
                    end
                case 'Custom (Subcohorts)'
                    if isempty(app.TrainOnListBox.Value) || isempty(app.PredictListBox.Value)
                        app.RunCVButton.Enable = 'off';
                    end
                    cvlEnable='on';
                otherwise
                    cvlEnable='off';
            end

            % CV Listboxes
            app.PredictListBox.Enable = cvlEnable;
            app.PredictListBoxLabel.Enable = cvlEnable;
            app.TrainOnListBox.Enable = cvlEnable;
            app.TrainOnListBoxLabel.Enable = cvlEnable;
            app.ModeltypeDropDown.Enable = cvlEnable;

            app.AutoRefreshCheckBox.Enable = val;
            app.RefreshViewButton.Enable = val;
            app.LiveVisualizeCheckBox.Enable = val;

            if ismember(app.Corrtype.Value,{'Skipped Spearman','Skipped Pearson'})
                app.ConsiderSignificantFibersonlyCheckBox.Enable = 'off';
                app.ConsiderSignificantFibersonlyCheckBox.Value = 0;
                app.levelEditField.Enable = 'off';
                app.levelEditFieldLabel.Enable = 'off';
                app.MultTestStrategy.Enable = 'off';
            else
                app.ConsiderSignificantFibersonlyCheckBox.Enable = val;
                app.levelEditField.Enable = val;
                app.levelEditFieldLabel.Enable = val;
                app.MultTestStrategy.Enable = val;
            end

            if ismember(app.Corrtype.Value,{'Skipped Spearman','Skipped Pearson','Bend'})
                app.CovariatesCheckBox.Enable = 'off';
                app.CovariatesCheckBox.Value = 0;
            end

            app.BasePredictiononDropDown.Enable = val;
            app.BasePredictiononDropDownLabel.Enable = val;
            %app.ModelSetupOptionButtonGroup.Enable = val;


            app.FittoscoresCheckBox.Enable = cvEnable;

            if app.SplitColorByPCAButton.Value % when doing a PCA based prediction, we need to fit back tract estimates to responsevars.
                app.FittoscoresCheckBox.Value = 1;
                app.FittoscoresCheckBox.Enable = 'off';
                app.tractset.doactualprediction = 1; % redundancy just to make sure in fact this should not be needed.
            end

            switch app.FittoscoresCheckBox.Value
                case 1
                    app.ModeltypeDropDown.Enable = cvEnable;
                case 0
                    app.ModeltypeDropDown.Enable = 'off';
            end

            if app.CovariatesCheckBox.Value
                app.CovariatesListBox.Enable = val;
            else
                app.CovariatesListBox.Enable = 'off';
            end

            switch app.ModelSetupDropDown.Value
                case {'T-Tests / VTAs (Baldermann 2019) / PAM (OSS-DBS)','T-Tests / PAM (OSS-DBS)'}
                    app.EfieldthresholdSlider.Enable = 'off';
                    app.EfieldthresholdLabel.Enable = 'off';
                    app.Efieldmetric.Enable = 'off';
                    app.Corrtype.Enable = 'off';
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'on';
                    app.BasePredictiononDropDown.Items={'Mean of Scores','Sum of Scores','Peak of Scores','Peak 5% of Scores','Histogram','z-scored Histogram'};
                    if ~ismember(app.tractset.basepredictionon,lower(app.BasePredictiononDropDown.Items))
                        app.tractset.basepredictionon=lower(app.BasePredictiononDropDown.Value);
                    end
                case 'Correlations / E-fields (Irmen 2020)'
                    app.EfieldthresholdSlider.Enable = val;
                    app.EfieldthresholdLabel.Enable = val;
                    app.Efieldmetric.Enable = val;
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'on';
                    app.Corrtype.Enable = val;
                    app.BasePredictiononDropDown.Items={'Mean of Scores','Sum of Scores','Peak of Scores','Peak 5% of Scores','Profile of Scores: Spearman','Profile of Scores: Pearson','Profile of Scores: Bend','Histogram','z-scored Histogram'};
                case 'Proportion Test (Chi-Square) / VTAs (binary vars)'
                    app.EfieldthresholdSlider.Enable = 'off';
                    app.EfieldthresholdLabel.Enable = 'off';
                    app.Efieldmetric.Enable = 'off';
                    app.Corrtype.Enable = 'off';
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'on';
                    app.BasePredictiononDropDown.Items={'Mean of Scores','Sum of Scores','Peak of Scores','Peak 5% of Scores','Histogram','z-scored Histogram'};
                    if ~ismember(app.tractset.basepredictionon,lower(app.BasePredictiononDropDown.Items))
                        app.tractset.basepredictionon=lower(app.BasePredictiononDropDown.Value);
                    end
                case 'Binomial Tests / VTAs (binary vars)'
                    app.EfieldthresholdSlider.Enable = 'off';
                    app.EfieldthresholdLabel.Enable = 'off';
                    app.Efieldmetric.Enable = 'off';
                    app.Corrtype.Enable = 'off';
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'on';
                    app.BasePredictiononDropDown.Items={'Mean of Scores','Sum of Scores','Peak of Scores','Peak 5% of Scores','Histogram','z-scored Histogram'};
                    if ~ismember(app.tractset.basepredictionon,lower(app.BasePredictiononDropDown.Items))
                        app.tractset.basepredictionon=lower(app.BasePredictiononDropDown.Value);
                    end
                case 'Reverse T-Tests / E-Fields (binary vars)'
                    app.EfieldthresholdSlider.Enable = val;
                    app.EfieldthresholdLabel.Enable = val;
                    app.Efieldmetric.Enable = val;
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'on';
                    app.Corrtype.Enable = val;
                    app.BasePredictiononDropDown.Items={'Mean of Scores','Sum of Scores','Peak of Scores','Peak 5% of Scores','Profile of Scores: Spearman','Profile of Scores: Pearson','Profile of Scores: Bend','Histogram','z-scored Histogram'};
                case 'Plain Connections'
                    app.EfieldthresholdSlider.Enable = 'off';
                    app.EfieldthresholdLabel.Enable = 'off';
                    app.Efieldmetric.Enable = 'off';
                    app.CovariatesCheckBox.Value = 0;
                    app.CovariatesCheckBox.Enable = 'off';
                    app.CovariatesListBox.Enable = 'off';
                    app.Corrtype.Enable = 'off';
                    app.ConsiderSignificantFibersonlyCheckBox.Value = 0;
                    app.ConsiderSignificantFibersonlyCheckBox.Enable = 'off';
                    app.levelEditField.Enable = 'off';
                    app.MultTestStrategy.Enable = 'off';
                    app.ShowNegativeFibersCheckBox.Value = 0;
                    app.ShowNegativeFibersCheckBox.Enable = 'off';
                    app.BasePredictiononDropDown.Items={'Sum of Scores','Mean of Scores','Peak of Scores','Peak 5% of Scores'};
                    if ~ismember(app.tractset.basepredictionon,lower(app.BasePredictiononDropDown.Items))
                        app.tractset.basepredictionon=lower(app.BasePredictiononDropDown.Value);
                    end
                    app.PosthoccorrectforgroupCheckBox.Value = 0;
            end


            switch app.StrategyDropDown.Value
                case 'Import Model'
                    app.tractset.useExternalModel = true;
                    app.ModelSetupDropDown.Enable = 'off';
                    app.ModelSetupDropDownLabel.Enable = 'off';
                    app.Efieldmetric.Enable = 'off';
                    app.Corrtype.Enable = 'off';

                    app.PredictListBoxLabel.Position(1) = 149;
                    app.PredictListBox.Position(1) = 200;
                    app.PredictListBox.Position(3) = 251;
                    app.TrainOnListBox.Visible = 'off';  % this might need to be placed above
                    app.ImportedModelName.Visible = 'on';

                otherwise
                    app.tractset.useExternalModel = false;
                    app.PredictListBoxLabel.Position(1) = 245;
                    app.PredictListBox.Position(1) = 297;
                    app.PredictListBox.Position(3) = 156;
                    app.TrainOnListBox.Visible = 'on';
                    app.ImportedModelName.Visible = 'off';
            end

            app.PatientsListBox.Enable=val;

            app.SubcohortsDropDown.Enable=val;
            app.SubcohortsDropDown.Items=[{'Choose...', 'All Patients', 'Create Subcohort from Selection', 'Create Subcohort from Inverse Selection'},app.tractset.setlabels];

            app.SubcohortsDropDown_2.Items=[app.tractset.setlabels];

            % Thresholding Method:
            app.ThresholdingDropDown.Enable = val;
            app.ThresholdingDropDownLabel.Enable = val;

            switch app.ShowNegativeFibersCheckBox.Value
                case 0
                    app.ShowAmountEditFieldNegative.Enable = 'off';
                    app.ColorButtonNegative.Enable = 'off';
                    app.ShowAmountLabelNegative.Enable = 'off';

                    app.ShowAmountLeftEditFieldNegative.Enable = 'off';
                    app.ShowAmountRightEditFieldNegative.Enable = 'off';

                    app.ShowAmountLeftLabelNegative.Enable = 'off';
                    app.ShowAmountRightLabelNegative.Enable = 'off';
                case 1
                    app.ShowAmountEditFieldNegative.Enable = val;
                    app.ColorButtonNegative.Enable = val;
                    app.ShowAmountLabelNegative.Enable = val;
                    if numsides==2 % usual bilateral case
                        app.ShowAmountLeftEditFieldNegative.Enable = val;
                        app.ShowAmountRightEditFieldNegative.Enable = val;

                        app.ShowAmountLeftLabelNegative.Enable = val;
                        app.ShowAmountRightLabelNegative.Enable = val;
                    else % scripting case e.g. when using only one side
                        app.ShowAmountLeftEditFieldNegative.Enable = 'off';
                        app.ShowAmountRightEditFieldNegative.Enable = 'off';

                        app.ShowAmountLeftLabelNegative.Enable = 'off';
                        app.ShowAmountRightLabelNegative.Enable = 'off';
                    end
            end

            switch app.ShowPositiveFibersCheckBox.Value
                case 0
                    app.ShowAmountEditFieldPositive.Enable = 'off';
                    app.ColorButtonPositive.Enable = 'off';
                    app.ShowAmountLabelPositive.Enable = 'off';

                    app.ShowAmountLeftEditFieldPositive.Enable = 'off';
                    app.ShowAmountRightEditFieldPositive.Enable = 'off';

                    app.ShowAmountLeftLabelPositive.Enable = 'off';
                    app.ShowAmountRightLabelPositive.Enable = 'off';
                case 1
                    app.ShowAmountEditFieldPositive.Enable = val;
                    app.ColorButtonPositive.Enable = val;
                    app.ShowAmountLabelPositive.Enable = val;
                    if numsides==2 % usual bilateral case
                        app.ShowAmountLeftEditFieldPositive.Enable = val;
                        app.ShowAmountRightEditFieldPositive.Enable = val;

                        app.ShowAmountLeftLabelPositive.Enable = val;
                        app.ShowAmountRightLabelPositive.Enable = val;
                    else % scripting case e.g. when using only one side
                        app.ShowAmountLeftEditFieldPositive.Enable = 'off';
                        app.ShowAmountRightEditFieldPositive.Enable = 'off';

                        app.ShowAmountLeftLabelPositive.Enable = 'off';
                        app.ShowAmountRightLabelPositive.Enable = 'off';
                    end

            end

            % Logic of multitract panel radio buttons - enabling of pca panel is
            % handled in updatepca()
            if app.SingleTractAnalysisButton.Value
                app.ColorButtonPositive.Enable = val;
                app.ColorButtonNegative.Enable = val;
                app.SubscoresListBox.Enable = 'off';
                app.CovariatesCheckBox.Enable = val;
                app.CovariatesListBox.Enable = val;
                app.MultitractWeightsButton.Enable = 'off';
                %app.SingleTractAnalysisTab.Parent = app.TabGroup;
                app.SingleTractAnalysisTab.Title = 'Single Tract Analysis';
                for i=1:10
                    tab_title = (['Subscore_' num2str(i) 'Tab']);
                    app.(tab_title).Parent = [];
                end
            end
            if app.SplitColorByGroupButton.Value
                app.ColorButtonPositive.Enable = 'off';
                app.ColorButtonNegative.Enable = 'off';
                app.SubscoresListBox.Enable = 'off';
                app.CovariatesCheckBox.Enable = val;
                app.CovariatesListBox.Enable = val;
                app.MultitractWeightsButton.Enable = 'off';
                n_groups = unique(app.tractset.M.patient.group);
                app.SingleTractAnalysisTab.Title = 'Group Fiber Visualization';
                for i=1:10
                    tab_title = (['Subscore_' num2str(i) 'Tab']);
                    if i <= length(n_groups)
                        app.(tab_title).Parent = app.TabGroup;
                        new_tabtitle = ['Group' num2str(i)];
                        if strcmp(app.(tab_title).Title,new_tabtitle)
                            app.(tab_title).Title = ['dummy' num2str(randi(2^52))];
                            app.(tab_title).Title = ['Group' num2str(i)];
                        else
                            app.(tab_title).Title = ['Group' num2str(i)];
                        end
                    else
                        app.(tab_title).Parent = [];
                    end
                end
            end
            if app.SplitColorBySubscoreButton.Value
                app.ColorButtonPositive.Enable = 'off';
                app.ColorButtonNegative.Enable = 'off';
                app.SubscoresListBox.Enable = val;
                app.CovariatesCheckBox.Enable = 'on';
                app.CovariatesListBox.Enable = 'on';
                app.MultitractWeightsButton.Enable = cvEnable;
                app.OpeninCleartuneButton.Enable = cvEnable;
                app.SingleTractAnalysisTab.Title = 'Mixed Fiber Visualization';
                for i=1:10
                    tab_title = (['Subscore_' num2str(i) 'Tab']);
                    if i <= length(app.tractset.subscore.labels)
                        app.(tab_title).Parent = app.TabGroup;
                        if strcmp(app.(tab_title).Title,app.tractset.subscore.labels{i})
                            %if you want to change the name to a name that
                            %is already present in MATLAB'S memory, it
                            %switches to the default name for this tab
                            %"Subscore_n"
                            app.(tab_title).Title = ['dummy' num2str(randi(2^52))];
                            app.(tab_title).Title = app.tractset.subscore.labels{i};
                        else
                            app.(tab_title).Title = app.tractset.subscore.labels{i};
                        end
                    else
                        app.(tab_title).Parent = [];
                    end
                end
            end
            if app.SplitColorByPCAButton.Value
                app.NestedCheckBox.Visible = 'off';
                app.ColorButtonPositive.Enable = val;
                app.ColorButtonNegative.Enable = val;
                app.SubscoresListBox.Enable = val;
                app.CovariatesCheckBox.Enable = 'off';
                app.CovariatesListBox.Enable = 'off';
                app.MultitractWeightsButton.Enable = 'off';
                app.SingleTractAnalysisTab.Title = 'PCA Fiber Visualization';
                for i=1:10
                    tab_title = (['Subscore_' num2str(i) 'Tab']);
                    if i <= app.tractset.numpcs
                        app.(tab_title).Parent = app.TabGroup;
                        new_tabtitle = ['PC' num2str(i)];
                        if strcmp(app.(tab_title).Title,new_tabtitle)
                            app.(tab_title).Title = ['dummy' num2str(randi(2^52))];
                            app.(tab_title).Title = ['PC' num2str(i)];
                        else
                            app.(tab_title).Title = ['PC' num2str(i)];
                        end
                    else
                        app.(tab_title).Parent = [];
                    end
                end
            end
            %end
            app.updatepca;
            % turn light back to green
            app.isbusy(0);
        end