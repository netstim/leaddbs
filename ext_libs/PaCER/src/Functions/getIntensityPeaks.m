%% Get peaks respectivly center points at threshold below the peaks from 1D
% intensity profile. Expected peakWaveCenter number == 4 for now (needs update to support electrodes with more ring contacts)!
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2014 - 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function [peakLocs, peakWaveCenters, peakValues, threshIntensityProfile, threshold, contactAreaCenter, contactAreaWidth, xrayMarkerAreaCenter, xrayMarkerAreaWidth] = getIntensityPeaks(intensityProfile, skelScaleMm, filterIdxs)
[peakValues, peakLocs, ~, pkPromineces] =  findpeaks(intensityProfile(filterIdxs),skelScaleMm(filterIdxs), 'MinPeakDistance', 1.4, 'MinPeakHeight', 1.1 * nanmean(intensityProfile), 'MinPeakProminence', 0.01 * nanmean(intensityProfile)); % find elec contacts, TODO 1.5 is for 3389/3387
xrayMarkerAreaWidth = [];
xrayMarkerAreaCenter = [];
try
    % find peakWaveCenters for each of the (four) expected single peaks
    threshold = min(peakValues(1:4)) - (min(pkPromineces(1:4)) / 4); % cool
    threshIntensityProfile = min(intensityProfile, threshold);
    contactSampleLabels = bwlabel(~(threshIntensityProfile(filterIdxs) < threshold));
    values = accumarray((contactSampleLabels+1)', skelScaleMm(filterIdxs)); % in x ([mm])
    counts = accumarray((contactSampleLabels+1)', 1);
    peakWaveCenters = values(2:5)./counts(2:5); % index 1 is the "zero label"
catch
   disp('peakWaveCenter detection failed. Returing peaksLocs in peakWaveCenters.');
   peakWaveCenters = peakLocs; 
   threshIntensityProfile = intensityProfile;
   threshold = NaN;
end

%% Detect "contact area" as fallback for very low SNR signals where no single contacts are visible
thresholdArea = mean(intensityProfile(filterIdxs));
threshIntensityProfileArea = min(intensityProfile, thresholdArea);
contactSampleLabels = bwlabel(~(threshIntensityProfileArea(filterIdxs) < thresholdArea));
values = accumarray((contactSampleLabels+1)', skelScaleMm(filterIdxs)); % in x ([mm])
counts = accumarray((contactSampleLabels+1)', 1);
contactAreaCenter = values(2)./counts(2); % index 1 is the "zero label", index 2  (value 1) is the contact region, index 3 (value 2) might be an X-Ray obaque arker
idxs = find(contactSampleLabels+1==2);
contactAreaWidth =  abs(skelScaleMm(idxs(1))-skelScaleMm(idxs(end)));
if(max(contactSampleLabels) > 1)
    disp('Multiple metal areas found along electrode. Is this an electrode type with an addtional X-Ray marker?');
    xrayMarkerAreaCenter = values(3)./counts(3); % index 1 is the "zero label", index 2  (value 1) is the contact region, index 3 (value 2) might be an X-Ray obaque arker
    idxs = find(contactSampleLabels+1==3);
    xrayMarkerAreaWidth =  abs(skelScaleMm(idxs(1))-skelScaleMm(idxs(end)));
end
%% Plot
% figure
% findpeaks(intensityProfile(filterIdxs),skelScaleMm(filterIdxs), 'MinPeakDistance', 3, 'MinPeakHeight', 1.1 * nanmean(intensityProfile), 'MinPeakProminence', 0.1 * nanmean(intensityProfile)); % find elec contacts, TODO 1.5 is for 3389/3387
% hold on;
% plot(skelScaleMm(filterIdxs), threshIntensityProfileArea(filterIdxs));
% scatter(peakWaveCenters, repmat(threshold,1,4), 'filled');
end
