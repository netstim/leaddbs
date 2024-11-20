function [vatlist,FilesExist] = ea_discfibers_getvats(obj)
% Return list of VATs


numPatient = length(obj.allpatients);
vatlist = cell(numPatient*2,2);
FilesExist = zeros(numPatient*2,2); % also check if mirrored fields exist

modelLabel = ea_simModel2Label(obj.M.vatmodel);

disp('Construct VAT list...')
for sub=1:numPatient

    [~, subj_tag] = fileparts(obj.allpatients{sub});

    % Original VAT E-field
    stimFolder = [obj.allpatients{sub}, filesep, 'stimulations', filesep, ea_nt(0), 'gs_', obj.M.guid];
    stimParams = [stimFolder, filesep, subj_tag, '_desc-stimparameters.mat'];

    if ~isfile(stimParams)
        FilesExist([sub,sub+numPatient], :) = ones(2);
        ea_cprintf('CmdWinWarnings', 'Stimulation parameters not found! Skip checking stimulation/vta existence.\n');
        return
    else
        S = ea_loadstimulation(stimParams);

        for side = 1:2
            % if no stimulation, we do not expect the file to exist
            % so no re-simulation is needed
            if all(S.amplitude{1,side} == 0)
                FilesExist(sub,side) = 1;
                vatlist{sub,side} = 'skip';

                % also mark mirrored
                if side == 1
                    vatlist{sub+numPatient,2} = 'skip';
                    FilesExist(sub+numPatient,2) = 1;
                else
                    vatlist{sub+numPatient,1} = 'skip';
                    FilesExist(sub+numPatient,1) = 1;
                end

                disp(subj_tag)
                disp(side)

                continue
            else
                if side == 1
                    BIDS_side = 'R';
                else
                    BIDS_side = 'L';
                end

                try
                    vatlist(sub,side) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-',BIDS_side,'.nii'], 0);
                    FilesExist(sub,side)=1;
                catch
                    if side == 1
                        ea_cprintf('CmdWinWarnings', 'Right side VTA doesn''t exist under stimulation folder:\n%s\n\n', stimFolder);
                    else
                        ea_cprintf('CmdWinWarnings', 'Left side VTA doesn''t exist under stimulation folder:\n%s\n\n', stimFolder);
                    end
                    vatlist(sub,side) = {''};
                    FilesExist(sub,side)=0;
                end
            end
        end
    end

    % Mirrored VAT E-field
    ea_genflippedjointnii(vatlist{sub,1}, vatlist{sub,2});
    try
        vatlist(numPatient+sub, 1) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-R_hemidesc-FlippedFromLeft\.nii$'], 0);
        FilesExist(numPatient+sub, 1)=1;
    catch
        if ~strcmp(vatlist(numPatient+sub, 2), "skip")
            ea_cprintf('CmdWinWarnings', 'Right side VTA (flipped from left) doesn''t exist under stimulation folder:\n%s\n\n', stimFolder);
            vatlist(numPatient+sub, 1) = {''};
        end
        %FilesExist(numPatient+sub, 1)=0;
    end
    try
        vatlist(numPatient+sub, 2) = ea_regexpdir(stimFolder, ['sim-efield_model-',modelLabel,'_hemi-L_hemidesc-FlippedFromRight\.nii$'], 0);
        FilesExist(numPatient+sub, 2)=1;
    catch
        if ~strcmp(vatlist(numPatient+sub, 2), "skip")
            ea_cprintf('CmdWinWarnings', 'Left side VTA (flipped from right) doesn''t exist under stimulation folder:\n%s\n\n', stimFolder);
            vatlist(numPatient+sub, 2) = {''};
        end
        %FilesExist(numPatient+sub, 2)=0;
    end
end

FilesExist=double(logical(FilesExist)); % convert to zeros and ones.
