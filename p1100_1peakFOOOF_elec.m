% 09/24/2024 Makoto. Modified.
% 05/30/2024 Makoto. Showing three results: 0, 1, and 3 peak models.
% 05/29/2024 Makoto. Fitting 3 Gaussians allowing broad speaks is not optimal for evaluating 1/f. 
% 01/10/2024 Makoto. Error detected.
% 01/09/2024 Makoto. Used.
% 01/05/2024 Makoto. Created.
clear
close all
clc

addpath(genpath('/srv/Makoto/Ribbon/code/ribbonanalysistoolbox'))

allSets          = dir('/srv/Makoto/TMS_depression/p0200_renameFiles/*.set');
allSetsProcessed = dir('/srv/Makoto/TMS_depression/p1000_concatenateDaySeparatedSessions/concatenated/*.set');

setNames = {allSets.name}';
subjIDs = cellfun(@(x) x(1:2), setNames, 'UniformOutput', false);
uniqueSubjIDs = unique(subjIDs);

goodFreqsHz = 2:0.1:52.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FOOOF parameter settting. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy', 'both'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       gaussian_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].


opt.freq_range          = [goodFreqsHz(1) goodFreqsHz(end)];
opt.aperiodic_mode      = 'fixed';   % aperiodic_mode : {'fixed', 'knee'} Defines absence or presence of knee in aperiodic component.
opt.max_peaks           = 1;        % Maximum number of gaussians to fit within the spectrum.
opt.peak_threshold      = 1;        % 2 std dev: parameter for interface simplification
opt.min_peak_height     = 0;        % Minimum height of a peak (in log10).
opt.peak_width_limits   = [0.5 3];    % Values were taken from https://neuroimage.usc.edu/brainstorm/Tutorials/Fooof
opt.proximity_threshold = 2;        % Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%opt.peak_type           = 'cauchy'; % {'gaussian', 'cauchy', 'both'}
opt.peak_type           = 'best'; % {'gaussian', 'cauchy', 'both'}
opt.guess_weight        = 'weak';   % Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
hOT                     = 1;        % Defines whether to use constrained optimization, fmincon, or basic simplex, fminsearch.
opt.thresh_after        = true;     % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)


idList = cell(1,1);
counter = 0;
for uniqueSubjIdx = 1:length(uniqueSubjIDs)

    currentID = uniqueSubjIDs{uniqueSubjIdx};

    currentIdMask = strcmp(subjIDs, currentID);

    baselineMask = contains(setNames, 'baseline');
    post1wkMask  = contains(setNames, 'post_1wk');
    post4wkMask  = contains(setNames, 'post_4wk');
    post12wkMask = contains(setNames, 'post_12wk');

    % Load baseline, 1wk, 4wk, 12wk datasets.
    baselineData = pop_loadset('filename', setNames(find(currentIdMask & baselineMask)), 'filepath', '/srv/Makoto/TMS_depression/p0200_renameFiles/');
    post1wkData  = pop_loadset('filename', setNames(find(currentIdMask & post1wkMask)),  'filepath', '/srv/Makoto/TMS_depression/p0200_renameFiles/');
    post4wkData  = pop_loadset('filename', setNames(find(currentIdMask & post4wkMask)),  'filepath', '/srv/Makoto/TMS_depression/p0200_renameFiles/');
    post12wkData = pop_loadset('filename', setNames(find(currentIdMask & post12wkMask)), 'filepath', '/srv/Makoto/TMS_depression/p0200_renameFiles/');

    if isempty(baselineData) | isempty(post1wkData) | isempty(post4wkData) | isempty(post12wkData)
        continue
    end

    % Load the concatenated process results.
    currentSubjBaselineSetName = setNames(find(currentIdMask & baselineMask));
    currentSubjIdx = find(contains({allSetsProcessed.name}, currentSubjBaselineSetName{1}(1:2)));
    EEGconcatenated = pop_loadset('filename', allSetsProcessed(currentSubjIdx).name, 'filepath', '/srv/Makoto/TMS_depression/p1000_concatenateDaySeparatedSessions/concatenated');
    
    EEG_cells{1,1} = baselineData;
    EEG_cells{2,1} = post1wkData;
    EEG_cells{3,1} = post4wkData;
    EEG_cells{4,1} = post12wkData;
    for conditionIdx = 1:4

        EEG = EEG_cells{conditionIdx};
        EEG.icaweights = EEGconcatenated.icaweights;
        EEG.icasphere  = EEGconcatenated.icasphere;
        EEG.icawinv    = EEGconcatenated.icawinv;
        EEG.icaact     = EEGconcatenated.icaweights*EEGconcatenated.icasphere*double(EEG.data);
        EEG.dipfit     = EEGconcatenated.dipfit;
        EEG.etc.ic_classification = EEGconcatenated.etc.ic_classification;

        %%%%%%%%%%%%%
        %%% FOOOF %%%
        %%%%%%%%%%%%%
        % Mean and Median PSD
        spectra_median = zeros(EEG.nbchan, 506);
        for chIdx = 1:EEG.nbchan
            currentChData = double(EEG.data(chIdx,:));

            % Run spectrogram() with the same parameters as above.
            [S,psdFreqs,T,P] = spectrogram(currentChData, EEG.srate, EEG.srate/2, 2:0.1:52.5, EEG.srate);
            psd_twoSided_mean   = 10*log10(mean(P,2)*2);
            psd_twoSided_median = 10*log10(median(P,2)*2);
            spectra_mean(  chIdx,:) = psd_twoSided_mean;
            spectra_median(chIdx,:) = psd_twoSided_median;
        end

        fooofOutputs_mean   = cell(size(spectra_mean,  1),1);
        fooofOutputs_median = cell(size(spectra_median,1),1);
        for chIdx = 1:EEG.nbchan

            %%%%%%%%%%%%%%%%%
            %%% Mean PSD. %%%
            %%%%%%%%%%%%%%%%%
            currentPSD = spectra_mean(chIdx,:);

            % Convert dB to microVolt^2 (the FOOOF requirement)
            currentPSD_microVoltSquare = 10.^(currentPSD/10);
            
            % Run FOOOF.
            [fs, fg] = FOOOF_matlab(currentPSD_microVoltSquare, goodFreqsHz, opt, hOT); % currentPSD MUST be a row vector.
            fg.fs    = fs;

            % Store FOOOF results.
            fooofOutputs_mean{chIdx} = fg;


            %%%%%%%%%%%%%%%%%%%
            %%% Median PSD. %%%
            %%%%%%%%%%%%%%%%%%%
            currentPSD = spectra_median(chIdx,:);

            % Convert dB to microVolt^2 (the FOOOF requirement)
            currentPSD_microVoltSquare = 10.^(currentPSD/10);
            
            % Run FOOOF.
            [fs, fg] = FOOOF_matlab(currentPSD_microVoltSquare, goodFreqsHz, opt, hOT); % currentPSD MUST be a row vector.
            fg.fs    = fs;

            % Store FOOOF results.
            fooofOutputs_median{chIdx} = fg;
        end
        EEG.etc.FOOOF_mean   = fooofOutputs_mean;
        EEG.etc.FOOOF_median = fooofOutputs_median;
        EEG.etc.PSD_concatenated.spectra_mean   = spectra_mean;
        EEG.etc.PSD_concatenated.spectra_median = spectra_median;
        EEG.etc.PSD_concatenated.freqs = psdFreqs;

        pop_saveset(EEG, 'filename', EEG.filename, 'filepath', '/srv/Makoto/TMS_depression/p1100_1peakFOOOF_elec')
    end
end