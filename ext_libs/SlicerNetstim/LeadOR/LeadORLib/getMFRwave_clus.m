function mfr = getMFRwave_clus(fileName, sr)

% fileName='C:\Users\Simon\Desktop\SBN6E96E\AO_2020-12-04_11-44-21_9.797000.csv';
% sr=44e3;

sig = csvread(fileName,1,0);

prevPath = pwd;
tmpPath = fullfile(fileparts(fileName), 'tmp');
mkdir(tmpPath);
cd(tmpPath);

mfr = [];

for i = 1:size(sig,2)
    
    tmpFile = ['for_wave_clus_' num2str(i) '.mat'];
    
    data = sig(:,i);
    
    save(tmpFile,'sr','data');
    
    Get_spikes({tmpFile});
    Do_clustering('all','make_plots',false);
    
    if ~isfile(['times_' tmpFile])
        load([tmpFile(1:end-4) '_spikes.mat'])
        mfr(i) = size(spikes,1)  / (length(data) / sr);
        continue
    end
    
    load(['times_' tmpFile])

    cluster_n = unique(cluster_class(:,1));
    cluster_n(cluster_n == 0) = [];
    cluster_spike_count = arrayfun(@(x) sum(cluster_class(:,1) == x), cluster_n);

    mfr(i) = median(cluster_spike_count) / (length(data) / sr);

end

cd(prevPath);
rmdir(tmpPath,'s');
