% function score = ea_practicelocalizations(user_data, ground_truth)
%     
%     
%     userReco = load(user_data);
%     correctReco = load(ground_truth);
% 
%     NumElectrodes = length(userReco.reco.props);
%     
%     for i = 1 : NumElectrodes  
%         userMarkers = userReco.reco.mni.markers(i);
%         groundMarkers = correctReco.reco.mni.markers(i);
%         
%         matrix_userMarkers = cell2mat(struct2cell(userMarkers));
%         matrix_groundMarkers = cell2mat(struct2cell(groundMarkers));
% 
%        %transform(i) = mldivide(matrix_groundMarkers, matrix_userMarkers);
%        for j = 1:4
%            for k = 1:3
%                 square_difference = sum((matrix_groundMarkers(j,k)-matrix_userMarkers(j,k))^2);
%            end
%        end
%        mean_square_difference = square_difference/16;
%     end
%     
%     
%     score = mean_square_difference;
% 
%     
% end

function score = ea_practicelocalizations(user_outputFolder, ground_truth)

    output_folder_dir = BIDSFetcher(user_outputFolder).datasetDir;
    subj_id = cell2mat(BIDSFetcher(user_outputFolder).subjId);
    user_reco_file = output_folder_dir + "/derivatives/leaddbs/sub-"+subj_id+"/reconstruction/sub-"+subj_id+"_desc-reconstruction.mat";
    
    % Hard coded for now, but will eventually grab this directory from a
    % folder within leaddbs
%     ground_truth_folder = '/Users/savirmadan/Documents/Localizations/OSF/LeadDBSTrainingDataset/OSFReconstructionFiles';
%     ground_truth_reco_file = ground_truth_folder + "sub-"+subj_id+"_desc-reconstruction.mat"
%     ground_truth = '/Users/savirmadan/Documents/Localizations/LeadDBS/Output/Patient0304Output/derivatives/leaddbs/sub-CbctDbs0304/reconstruction/sub-CbctDbs0304_desc-reconstruction.mat';
    
    userReco = load(user_reco_file);
    
    correctReco = load(ground_truth);
%     correctReco = load(ground_truth_reco_file);
    
    score_array = [];
    square_difference = [];
    NumElectrodes = length(userReco.reco.props);
    
    for i = 1 : NumElectrodes  
        userMarkers = userReco.reco.mni.markers(i);
        groundMarkers = correctReco.reco.mni.markers(i);
        
        matrix_userMarkers = cell2mat(struct2cell(userMarkers));
        matrix_groundMarkers = cell2mat(struct2cell(groundMarkers));
        user_x_head_diff = matrix_userMarkers(3,:) - matrix_userMarkers(1,:);
        normalizedX = matrix_groundMarkers(1,:) + user_x_head_diff;
        matrix_userMarkers(3,:) = normalizedX;
       %transform(i) = mldivide(matrix_groundMarkers, matrix_userMarkers);
       for j = 1:4
           for k = 1:3
                square_difference(j) = sum((matrix_groundMarkers(j,k)-matrix_userMarkers(j,k))^2);
           end
           temp_square_difference = square_difference;
           square_difference(j) = sqrt(square_difference(j));
       end
       elecdiff{i} = square_difference;
       elec_score = matrix_groundMarkers - matrix_userMarkers;
       score_array{i} = elec_score;
    end
    
    % Temporary score
    score = elecdiff;
    fileID = fopen('score_output.txt','w');
    fprintf(fileID,'This file contains information on your electrode localization accuracy:\n');
    fprintf(fileID,'\nElectrode Marker Differences:\n');
    for i = 1:numel(elecdiff)
        fprintf(fileID, 'Electrode %d:\n', i);
%         fprintf(fileID, '%f ', elecdiff{i});
        for j = 1:3 
            if j == 1
                fprintf(fileID, 'Difference in head marker (mm): ');
            end
            if j == 2
                fprintf(fileID, 'Difference in tail marker (mm): ');
            end
            if j == 3
                fprintf(fileID, 'Difference in rotation marker (mm): ');

            end
            fprintf(fileID, '%f', elecdiff{i}(j));    
            fprintf(fileID, '\n');

        end
        fprintf(fileID, '\n');  

    end
    fclose(fileID); % Close the file after writing
    
    % Open the text file using the default application
    open('score_output.txt');
    x = ["head" "tail" "x"];
    y1 = [elecdiff{1}(1) elecdiff{1}(2) elecdiff{1}(3)];
    y2 = [elecdiff{2}(1) elecdiff{2}(2) elecdiff{2}(3)];
    Y = [y1
        y2];
    figure
    bar(Y);


end

    




