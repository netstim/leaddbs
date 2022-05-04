function [Ihat,Ihattrain] = compute_fibscore_model(obj, fibsval, Ihat, Ihattrain, patientsel, training, test, Iperm_dummy,Improvement,test_index)

    % args: vals, obj, I_hat, training, test, fibsval,
    % Iperm, patientsel, opt: Slope, Intercept
    corrtype = 'spearman';
    if ~exist('Iperm', 'var')
        if obj.cvlivevisualize
            [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
            obj.draw(vals,fibcell,usedidx)
            %obj.draw(vals,fibcell);
            drawnow;
        else
            [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training));
        end
    else
        if obj.cvlivevisualize
            [vals,fibcell,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
            obj.draw(vals,fibcell,usedidx)
            %obj.draw(vals,fibcell);
            drawnow;
        else
            [vals,~,usedidx] = ea_discfibers_calcstats(obj, patientsel(training), Iperm);
        end
    end

    if size(obj.responsevar,2)>1
        lateral_score = true;
    else
        lateral_score = false;
    end
    
    switch obj.modelNormalization
        case 'z-score'
            for s=1:length(vals)
                vals{s}=ea_nanzscore(vals{s});
            end
        case 'van Albada 2007'
            for s=1:length(vals)
                vals{s}=ea_normal(vals{s});
            end
    end
    

    for voter=1:size(vals,1)
        for side=1:size(vals,2)

            % if not lateral, compute Ihat for both hemispheres at the same
            % time
            if lateral_score == false                 
                fibsval_usedidx = cell(1,2);
                fibsval_usedidx{1} = fibsval{1,1}(usedidx{voter,1},patientsel(test));
                fibsval_usedidx{2} = fibsval{1,2}(usedidx{voter,2},patientsel(test));
                vals_flat = vertcat(vals{voter,:});
                fibsval_usedidx_test_flat = vertcat(fibsval_usedidx{:});
                
                fibsval_usedidx_train = cell(1,2);
                fibsval_usedidx_train{1} = fibsval{1,1}(usedidx{voter,1},patientsel(training));
                fibsval_usedidx_train{2} = fibsval{1,2}(usedidx{voter,2},patientsel(training));
                vals_flat_train = vertcat(vals{voter,:});
                fibsval_usedidx_train_flat = vertcat(fibsval_usedidx_train{:});
            end    

            if ~isempty(vals{voter,side})

                switch obj.statmetric % also differentiate between methods in the prediction part.
                    case {1,3,4} % ttests / OSS-DBS / reverse t-tests
                        switch lower(obj.basepredictionon)
                            case 'mean of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nanmean(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    Ihattrain(training,1,voter) = ea_nanmean(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihattrain(training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);                              
                                end

                            case 'sum of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nansum(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihattrain(training,1,voter) = ea_nansum(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihattrain(training,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                end
                            case 'peak of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nanmax(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    Ihattrain(training,1,voter) = ea_nanmax(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    Ihattrain(training,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                end
                            case 'peak 5% of scores'
                                if lateral_score == false
                                    ihatvals = vals_flat.*fibsval_usedidx_test_flat;
                                    ihatvals=sort(ihatvals);
                                    Ihat(test,1, voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 

                                    
                                    ihatvals_train = vals_flat_train.*fibsval_usedidx_train_flat;
                                    ihatvals_train=sort(ihatvals_train);
                                    Ihattrain(training,1,voter) = ea_nansum(ihatvals_train(1:ceil(size(ihatvals_train,1).*0.05),:),1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    break % both sides are already filled out!
                                else                                
                                    ihatvals=vals{1,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                    ihatvals=sort(ihatvals);
                                    Ihat(test,side,voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                    ihatvals_train=vals{1,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                    ihatvals_train=sort(ihatvals_train);
                                    Ihattrain(training,side,voter) = ea_nansum(ihatvals_train(1:ceil(size(ihatvals_train,1).*0.05),:),1);
                                end
                                
                                
                        end
                    case {2,5} % efields
                        switch lower(obj.basepredictionon)
                            case 'profile of scores: spearman'
                                %%%%TODO: not sure how train values are set
                                if lateral_score == false
                                    Ihat(test,1, voter) = corr(vals_flat,fibsval_usedidx_test_flat,'rows','pairwise','type','spearman');
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = (corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','spearman'));
                                end

                                if any(isnan(Ihat(test,side,voter)))
                                    Ihat(isnan(Ihat(test,side,voter)),side)=0;
                                    warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                end
                                Ihattrain(training,side,voter) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','spearman'));
                                
                            case 'profile of scores: pearson'
                                if lateral_score == false
                                    Ihat(test,1, voter) = corr(vals_flat,fibsval_usedidx_test_flat,'rows','pairwise','type','pearson');
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = (corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test)),'rows','pairwise','type','pearson'));
                                end
                                if any(isnan(Ihat(test,side,voter)))
                                    Ihat(isnan(Ihat(test,side,voter)),side)=0;
                                    warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                end
                                Ihattrain(training,side,voter) = atanh(corr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training)),'rows','pairwise','type','pearson'));
                                
                            case 'profile of scores: bend'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_bendcorr(vals_flat,fibsval_usedidx_test_flat);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = (ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(test))));
                                end
                                if any(isnan(Ihat(test,side,voter)))
                                    Ihat(isnan(Ihat(test,side,voter)),side)=0;
                                    warning('Profiles of scores could not be evaluated for some patients. Displaying these points as zero entries. Lower threshold or carefully check results.');
                                end
                                Ihattrain(training,side,voter) = atanh(ea_bendcorr(vals{voter,side},fibsval{1,side}(usedidx{voter,side},patientsel(training))));
                                
                            case 'mean of scores'
                                if ~isempty(vals{voter,side})
                                    if lateral_score == false
                                        Ihat(test,1, voter) = ea_nanmean(vals_flat.*fibsval_usedidx_test_flat,1);
                                        Ihat(test,2, voter) = Ihat(test,1, voter); 
                                        Ihattrain(training,1,voter) = ea_nanmean(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                        Ihattrain(training,2, voter) = Ihattrain(training,1, voter);
                                        break % both sides are already filled out!
                                    else          
                                        Ihat(test,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    end
                                end
                                Ihattrain(training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                
                            case 'sum of scores'
                                if ~isempty(vals{voter,side})
                                    if lateral_score == false
                                        Ihat(test,1, voter) = ea_nansum(vals_flat.*fibsval_usedidx_test_flat,1);
                                        Ihat(test,2, voter) = Ihat(test,1, voter); 
                                        Ihattrain(training,1,voter) = ea_nansum(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                        Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                        break % both sides are already filled out!
                                    else    
                                        Ihat(test,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    end
                                end
                                Ihattrain(training,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                
                            case 'peak of scores'
                                if ~isempty(vals{voter,side})
                                    if lateral_score == false
                                        Ihat(test,1, voter) = ea_nanmax(vals_flat.*fibsval_usedidx_test_flat,1);
                                        Ihat(test,2, voter) = Ihat(test,1, voter); 
                                        Ihattrain(training,1,voter) = ea_nanmax(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                        Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                        break % both sides are already filled out!
                                    else    
                                        Ihat(test,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                    end
                                end
                                Ihattrain(training,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                
                            case 'peak 5% of scores'
                                if ~isempty(vals{voter,side})
                                    if lateral_score == false
                                        ihatvals = vals_flat.*fibsval_usedidx_test_flat;
                                        ihatvals=sort(ihatvals); 
                                        ihatvals_train = vals_flat_train.*fibsval_usedidx_train_flat;
                                        ihatvals_train = sort(ihatvals_train);
                                        Ihat(test,side,voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                        Ihattrain(test,side,voter) = ea_nansum(ihatvals_train(1:ceil(size(ihatvals_train,1).*0.05),:),1);
                                        Ihat(test,2, voter) = Ihat(test,1, voter); 
                                        Ihattrain(training,2,voter) = Ihattrain(training,1,voter);
                                        break
                                    else    
                                        ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                        ihatvals=sort(ihatvals); 
                                        ihatvals_train=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                        ihatvals_train = sort(ihatvals_train);
                                        Ihattrain(training,side,voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals_train,1).*0.05),:),1);
                                        
                                    end
                                end
                                
                                
                        end

                    case 6 % Plain Connection
                        switch lower(obj.basepredictionon)
                            case 'mean of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nanmean(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    Ihattrain(training,1,voter) = ea_nanmean(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter);
                                    break % both sides are already filled out!
                                else                                
                                    Ihat(training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                end
                                Ihattrain(training,side,voter) = ea_nanmean(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                
                            case 'sum of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nansum(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter);
                                    Ihattrain(training,1,voter) = ea_nansum(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    break % both sides are already filled out!
                                else
                                    Ihat(test,side) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                end
                                Ihattrain(training,side,voter) = ea_nansum(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                                
                            case 'peak of scores'
                                if lateral_score == false
                                    Ihat(test,1, voter) = ea_nanmax(vals_flat.*fibsval_usedidx_test_flat,1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    Ihattrain(training,1,voter) = ea_nanmax(vals_flat_train.*fibsval_usedidx_train_flat,1);
                                    Ihattrain(training,2, voter) = Ihattrain(training,1, voter); 
                                    %break % both sides are already filled out!
                                else
                                    Ihat(test,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test)),1);
                                end
                                Ihattrain(training,side,voter) = ea_nanmax(vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training)),1);
                            case 'peak 5% of scores' %this peak 5%of scores is not consistent with the rest of the methods????
                                if lateral_score == false
                                    ihatvals = vals_flat.*fibsval_usedidx_test_flat;
                                    ihatvals=sort(ihatvals);
                                    Ihat(test,1, voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                    Ihat(test,2, voter) = Ihat(test,1, voter); 
                                    break % both sides are already filled out!
                                else  
                                    ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(test));
                                    ihatvals=sort(ihatvals);
                                    Ihat(test,side,voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                                end
                                ihatvals=vals{voter,side}.*fibsval{1,side}(usedidx{voter,side},patientsel(training));
                                ihatvals=sort(ihatvals);
                                Ihattrain(training,side,voter) = ea_nansum(ihatvals(1:ceil(size(ihatvals,1).*0.05),:),1);
                        end
                end
            end
        end
    end
%     Ihat = Ihat(test_index,:,:);
%     selected_pts = test;
%     weightmatrix=zeros(size(Ihat));
%     for voter=1:size(Ihat,3)
%         if ~isnan(obj.subscore.weights(voter)) % same weight for all subjects in that voter (slider was used)
%             weightmatrix(:,:,voter)=obj.subscore.weights(voter);
%         else % if the weight value is nan, this means we will need to derive a weight from the variable of choice
%             weightmatrix(:,:,voter)=repmat(obj.subscore.weightvars{voter}(selected_pts),1,size(weightmatrix,2)/size(obj.subscore.weightvars{voter}(selected_pts),2));
%         end
%     end
%     for xx=1:length(test_index)
%         for yy=1:size(Ihat,2)
%             %for xx=1:size(Ihat_voters,1) % make sure voter weights sum up to 1
%             %    for yy=1:size(Ihat_voters,2)
%             weightmatrix(xx,yy,:)=weightmatrix(xx,yy,:)./ea_nansum(weightmatrix(xx,yy,:));
%         end
%     end
%     Ihat=ea_nansum(Ihat.*weightmatrix,3);
%     Ihat = ea_nanmean(Ihat,2);
%     I = Improvement(test_index);
%     [R,p]=corr(I,Ihat,'rows','pairwise','type',corrtype);
%     fprintf("R value for set is %.2f\n",R);
%     fprintf("p value for set is %.2f\n",p);