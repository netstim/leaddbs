%% plotIntensityProfileAndPeaks
%
% Andreas Husch
% Centre Hospitalier de Luxembourg, Dep. of Neurosurgery /
% University of Luxembourg - Luxembourg Centre for Systems Biomedicne
% 2017
% mail@andreashusch.de, husch.andreas@chl.lu
function plotIntensityProfileAndPeaks(intensityProfile, skelScaleMm)
filterIdxs = find(skelScaleMm <= 20); % restrict to first 20mm
[peakLocs, peakWaveCenters, peakValues, threshIntensityProfile, threshold] = getIntensityPeaks(intensityProfile, skelScaleMm, filterIdxs);

figure('Name', 'Intensity Profile');
ax = gca;
xlabel('Length [mm]');
ylabel('Intensity Measure [HU]');
ylim([0 max(intensityProfile(filterIdxs))*1.05]);
xlim([0 20]) % show only first 20 mm
hold on;
hProfile = plot(skelScaleMm, intensityProfile);
% Add the threshold used to detect the ContactArea as a line
plot(skelScaleMm, ones(size(skelScaleMm))*mean(intensityProfile(filterIdxs)))
plot(peakLocs, peakValues + 0.03 * peakValues, 'v', 'MarkerFaceColor', hProfile.Color, 'MarkerEdgeColor', hProfile.Color);
ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
plot(skelScaleMm(filterIdxs), threshIntensityProfile(filterIdxs));
scatter(peakWaveCenters, repmat(threshold,1,length(peakWaveCenters)), 'filled');
grid on;
end