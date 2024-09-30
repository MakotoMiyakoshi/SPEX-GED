% 09/24/2024 Makoto. Created.
clear
close all
clc

addpath(genpath('/srv/Makoto/Ribbon/code/ribbonanalysistoolbox'))

allSets = dir('/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec/*.set');
setNames = {allSets.name}';
subjIDs = cellfun(@(x) x(1:2), setNames, 'UniformOutput', false);
uniqueSubjIDs = unique(subjIDs);

baselineMeanMatrix   = zeros(128, 8);
baselineMedianMatrix = zeros(128, 8);
post1wkMeanMatrix    = zeros(128, 8);
post1wkMedianMatrix  = zeros(128, 8);
post4wkMeanMatrix    = zeros(128, 8);
post4wkMedianMatrix  = zeros(128, 8);
post12wkMeanMatrix   = zeros(128, 8);
post12wkMeanMatrix   = zeros(128, 8);
for uniqueSubjIdx = 1:length(uniqueSubjIDs)

    currentID = uniqueSubjIDs{uniqueSubjIdx};

    currentIdMask = strcmp(subjIDs, currentID);

    baselineMask = contains(setNames, 'baseline');
    post1wkMask  = contains(setNames, 'post_1wk');
    post4wkMask  = contains(setNames, 'post_4wk');
    post12wkMask = contains(setNames, 'post_12wk');

    % Load baseline, 1wk, 4wk, 12wk datasets.
    baselineEEG = pop_loadset('filename', setNames(find(currentIdMask & baselineMask)), 'filepath', '/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec');
    post1wkEEG  = pop_loadset('filename', setNames(find(currentIdMask & post1wkMask)),  'filepath', '/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec');
    post4wkEEG  = pop_loadset('filename', setNames(find(currentIdMask & post4wkMask)),  'filepath', '/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec');
    post12wkEEG = pop_loadset('filename', setNames(find(currentIdMask & post12wkMask)), 'filepath', '/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec');

    spexVec_baseline_mean   = zeros(baselineEEG.nbchan,1);
    spexVec_baseline_median = zeros(baselineEEG.nbchan,1);
    for chIdx = 1:baselineEEG.nbchan
        currentSpex_mean   = baselineEEG.etc.FOOOF_mean{  chIdx}.aperiodic_params(2);
        currentSpex_median = baselineEEG.etc.FOOOF_median{chIdx}.aperiodic_params(2);
        spexVec_baseline_mean(chIdx)   = currentSpex_mean;
        spexVec_baseline_median(chIdx) = currentSpex_median;
    end
    baselineMeanMatrix(  :,uniqueSubjIdx) = spexVec_baseline_mean;
    baselineMedianMatrix(:,uniqueSubjIdx) = spexVec_baseline_median;
        %{
        figure
        plot(spexVec_baseline_mean)
        hold on
        plot(spexVec_baseline_median)
        %}

    spexVec_post1wk_mean   = zeros(baselineEEG.nbchan,1);
    spexVec_post1wk_median = zeros(baselineEEG.nbchan,1);
    for chIdx = 1:baselineEEG.nbchan
        currentSpex_mean   = post1wkEEG.etc.FOOOF_mean{  chIdx}.aperiodic_params(2);
        currentSpex_median = post1wkEEG.etc.FOOOF_median{chIdx}.aperiodic_params(2);
        spexVec_post1wk_mean(chIdx)   = currentSpex_mean;
        spexVec_post1wk_median(chIdx) = currentSpex_median;
    end
    post1wkMeanMatrix(  :,uniqueSubjIdx) = spexVec_post1wk_mean;
    post1wkMedianMatrix(:,uniqueSubjIdx) = spexVec_post1wk_median;

    spexVec_post4wk_mean   = zeros(baselineEEG.nbchan,1);
    spexVec_post4wk_median = zeros(baselineEEG.nbchan,1);
    for chIdx = 1:baselineEEG.nbchan
        currentSpex_mean   = post4wkEEG.etc.FOOOF_mean{  chIdx}.aperiodic_params(2);
        currentSpex_median = post4wkEEG.etc.FOOOF_median{chIdx}.aperiodic_params(2);
        spexVec_post4wk_mean(chIdx)   = currentSpex_mean;
        spexVec_post4wk_median(chIdx) = currentSpex_median;
    end
    post4wkMeanMatrix(  :,uniqueSubjIdx) = spexVec_post4wk_mean;
    post4wkMedianMatrix(:,uniqueSubjIdx) = spexVec_post4wk_median;

    spexVec_post12wk_mean   = zeros(baselineEEG.nbchan,1);
    spexVec_post12wk_median = zeros(baselineEEG.nbchan,1);
    for chIdx = 1:baselineEEG.nbchan
        currentSpex_mean   = post12wkEEG.etc.FOOOF_mean{  chIdx}.aperiodic_params(2);
        currentSpex_median = post12wkEEG.etc.FOOOF_median{chIdx}.aperiodic_params(2);
        spexVec_post12wk_mean(chIdx)   = currentSpex_mean;
        spexVec_post12wk_median(chIdx) = currentSpex_median;
    end
    post12wkMeanMatrix(  :,uniqueSubjIdx) = spexVec_post12wk_mean;
    post12wkMedianMatrix(:,uniqueSubjIdx) = spexVec_post12wk_median;
end


%% Plot single-subject SPEX topo for (mean vs. median) vs. (Base vs. 12 wk vs. Diff)
figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(baselineMeanMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, Base', subjIdx));
        colorbar;
    else
        topoplot(mean(baselineMeanMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/baseline_mean.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/baseline_mean', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(baselineMedianMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d median, Base', subjIdx));
        colorbar;
    else
        topoplot(mean(baselineMedianMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/baseline_median.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/baseline_median', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post12wkMeanMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 12wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post12wkMeanMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/post12wk_mean.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/post12wk_mean', '-djpeg95', '-r200')



figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post12wkMedianMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 12wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post12wkMedianMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/post12wk_median.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/post12wk_median', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post12wkMeanMatrix(:,subjIdx)-baselineMeanMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 12wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post12wkMeanMatrix-baselineMeanMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk12.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk12', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post12wkMedianMatrix(:,subjIdx)-baselineMedianMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 12wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post12wkMedianMatrix-baselineMedianMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_median_wk12.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_median_wk12', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post1wkMeanMatrix(:,subjIdx)-baselineMeanMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 1wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post1wkMeanMatrix-baselineMeanMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk1.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk1', '-djpeg95', '-r200')


figure('position', [50 50 1000 800])
tiledlayout(3,3,"TileSpacing", "tight")
for subjIdx = 1:9
    nexttile

    if subjIdx<=8
        topoplot(post4wkMeanMatrix(:,subjIdx)-baselineMeanMatrix(:,subjIdx), baselineEEG.chanlocs);
        title(sprintf('S%d mean, 4wk', subjIdx));
        colorbar;
    else
        topoplot(mean(post4wkMeanMatrix-baselineMeanMatrix,2), baselineEEG.chanlocs);
        title('Grand mean');
        colorbar;
    end
end
exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk4.pdf', 'ContentType', 'vector')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/diff_mean_wk4', '-djpeg95', '-r200')


%% GED applilcation.

covS = post12wkMeanMatrix*post12wkMeanMatrix';
covR = baselineMeanMatrix*baselineMeanMatrix';

% Generalized eigendecomposition (GED)
[evecs,evals] = eig(covS,covR);

% sort eigenvalues/vectors
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);

% plot the eigenspectrum
figure('position', [50 50 2100 800])
subplot(1,3,1)
plot(evals,'ks-','markersize',10,'markerfacecolor','m')
axis square
title('GED eigenvalues')
xlabel('Component number'), ylabel('Power ratio (\lambda)')

% The filter forward model is what the source "sees" when it looks through the
% electrodes. It is obtained by passing the covariance matrix through the filter.
subplot(1,3,2)
filt_topo = covS*evecs(:,1);
[~,se] = max(abs( filt_topo ));
filt_topo = filt_topo * sign(filt_topo(se));
filt_topo = filt_topo-mean(filt_topo);
topoplot(filt_topo, baselineEEG.chanlocs)
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '\DeltaSPEX')
title('GED Comp.1 (mean subtracted)')

subplot(1,3,3)
filt_topo = covS*evecs(:,2);
[~,se] = max(abs( filt_topo ));
filt_topo = filt_topo * sign(filt_topo(se));
filt_topo = filt_topo-mean(filt_topo);
filt_topo = filt_topo*-1; % Manual correction.
topoplot(filt_topo, baselineEEG.chanlocs)
cbarHandle = colorbar;
set(get(cbarHandle, 'title'), 'string', '\DeltaSPEX')
title('GED Comp.2 (mean subtracted)')

print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/comparisonForBamAward', '-djpeg95', '-r200')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/comparisonForBamAward', '-dsvg')
print('/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/comparisonForBamAward', '-dsvg')

exportgraphics(gcf, '/srv/Makoto/TMS_depression/p1110_GED_baseline_12weeks/comparisonForBamAward.pdf', 'ContentType', 'vector')
