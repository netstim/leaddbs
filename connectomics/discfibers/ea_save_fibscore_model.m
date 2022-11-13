function ea_save_fibscore_model(obj, patientsel, training, Iperm)

    % get fiber model (vals) and corresponding indices (usedidx)
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

    % maps fiber model (vals) into the full connectome space
    vals_all = cell(size(vals)); 
    for voter=1:size(vals,1)
        for side=1:size(vals,2)
            vals_all{voter,side} = zeros(obj.results.(ea_conn2connid(obj.connectome)).totalFibers,1);
            switch obj.connectivity_type
                case 2
                    vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_PAM{side}(usedidx{voter,side})) = vals{voter,side};
                otherwise
                    vals_all{voter,side}(obj.results.(ea_conn2connid(obj.connectome)).connFiberInd_VAT{side}(usedidx{voter,side})) = vals{voter,side};
            end
        end
    end

    % add the connectome name for recognition
    ftr.vals_all = vals_all;
    ftr.connectome = ea_conn2connid(obj.connectome);
    ftr.fibsvalType = ea_method2methodid(obj);
    save([obj.M.root,'vals_all_model.mat'], '-struct', 'ftr');

end
