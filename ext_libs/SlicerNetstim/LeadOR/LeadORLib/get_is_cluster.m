function is_cluster = get_is_cluster(recording_file)

recording_file = char(recording_file);

%%
is_cluster = false;

[file_path,file_name,~] = fileparts(recording_file);

tmp_path = fullfile(file_path,'tmp');
mkdir(tmp_path);
cd(tmp_path);

file_name = strrep(file_name,'.','_');

file_raw = fullfile(tmp_path, [file_name '.h5']);
file_spike = fullfile(tmp_path, [file_name '_spikes.mat']);
file_cluster = fullfile(tmp_path, ['times_' file_name '.mat']);
file_single_unit = fullfile(tmp_path, ['single_unit_' file_name '.mat']);

copyfile(recording_file, file_raw);

Get_spikes(cellstr({file_raw}), 'par', struct('detection','neg'));
if ~isfile(file_spike)
    finally(tmp_path)
    return
end

Do_clustering(cellstr({file_spike}), 'make_plots', false);
if ~isfile(file_cluster)
    finally(tmp_path)
    return
end

create_single_unit(file_cluster)
if ~isfile(file_single_unit)
    finally(tmp_path)
    return
end

wave_clus_load = load(file_single_unit);
[p2p, snr] = get_p2p_snr(wave_clus_load.spikes);
proportion_3ms = sum(diff(wave_clus_load.cluster_class(:,2)) < 3e-3) / size(wave_clus_load.cluster_class,1);

if (snr > 1.5)  && (p2p < 1e3)  &&  (proportion_3ms < 0.1) && (size(wave_clus_load.cluster_class,1) > 100)
    is_cluster = true;
end

finally(tmp_path)

end

function finally(tmp_path)
cd(fileparts(fileparts(tmp_path)));
rmdir(tmp_path,'s');
end

%%
function [] = create_single_unit(file_cluster)

load(file_cluster, 'cluster_class');
cluster_class(cluster_class(:,1)==0,:) = []; % remove 0 index
unique_cluster_index = unique(cluster_class(:,1));
unique_cluster_index = reshape(unique_cluster_index,1,[]);
    
snr = [];

for c = unique_cluster_index

    wave_clus_load = load(file_cluster);
    % keep data from this cluster
    idx = wave_clus_load.cluster_class(:,1) == c;
    wave_clus_load.cluster_class = wave_clus_load.cluster_class(idx,:);
    wave_clus_load.inspk         = wave_clus_load.inspk(idx,:);
    wave_clus_load.spikes        = wave_clus_load.spikes(idx,:);
    % set cluster idx to 1
    wave_clus_load.cluster_class(:,1) = 1;
    
    [~, snr(end+1)] = get_p2p_snr(wave_clus_load.spikes);
    
end

unique_cluster_index(snr<1.5) = [];

if isempty(unique_cluster_index)
    return
end

wave_clus_load = load(file_cluster);
% keep data from snr > 1.5
idx = any(wave_clus_load.cluster_class(:,1) == unique_cluster_index, 2);
wave_clus_load.cluster_class = wave_clus_load.cluster_class(idx,:);
wave_clus_load.inspk         = wave_clus_load.inspk(idx,:);
wave_clus_load.spikes        = wave_clus_load.spikes(idx,:);
% set cluster idx to 1
wave_clus_load.cluster_class(:,1) = 1;

out_fname = strrep(file_cluster, 'times', 'single_unit');
save(out_fname, '-struct', 'wave_clus_load')

end


function [p2p,snr] = get_p2p_snr(spikes)
mean_waveform = mean(spikes, 1);
p2p = max(mean_waveform) - min(mean_waveform);
snr = p2p / (std(reshape(spikes-mean_waveform, 1, [])) * 5);      
end

