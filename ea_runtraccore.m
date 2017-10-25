function [coords_mm,trajectory,markers]=ea_runtraccore(options)
directory = [options.root,options.patientname,filesep];


                    for side=options.sides

                        %try
                        % call main routine reconstructing trajectory for one side.
                        [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(options.patientname,options,side);

                        % refit electrodes starting from first electrode (this is redundant at this point).
                        coords_mm{side} = ea_map_coords(coords', [directory,options.prefs.tranii])';

                        [~,distmm]=ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),[directory,options.prefs.tranii]);

                        comp = ea_map_coords([0,0,0;trajvector{side}]', [directory,options.prefs.tranii])'; % (XYZ_mm unaltered)

                        trajvector{side}=diff(comp);

                        normtrajvector{side}=trajvector{side}./norm(trajvector{side});

                        for electrode=2:4
                            coords_mm{side}(electrode,:)=coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
                        end
                        markers(side).head=coords_mm{side}(1,:);
                        markers(side).tail=coords_mm{side}(4,:);

                        orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);

                        markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
                        markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality

                        coords_mm=ea_resolvecoords(markers,options);
                    end

                    % transform trajectory to mm space:
                    for side=1:length(options.sides)
                        try
                            if ~isempty(trajectory{side})
                                trajectory{side}=ea_map_coords(trajectory{side}', [directory,options.prefs.tranii])';
                            end

                        end
                    end
                    options.hybridsave=1;

                    % save reconstruction results
                    ea_methods(options,...
                        ['DBS-Electrodes were automatically pre-localized in native & template space using Lead-DBS software',...
                        ' (Horn & Kuehn 2015; SCR_002915; http://www.lead-dbs.org).'],...
                        {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});
             