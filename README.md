# SPEX-GED
Spectral exponent analysis with generalized eigenvalue decomposition

(Paper in prep)

Title: SPEX-GED: A whole-brain E-I balance estimation from EEG and contrast maximization

Authors: Hyeonseok Kim1, Makoto Miyakoshi2

Affiliations:
1 Swartz Center for Computational Neuroscience, University of California San Diego
2 Cincinnati Children’s Hospital Medical Center, University of Cincinnati College of Medicine

Abstract

Introduction
The goal of the SPEX-GED algorithm is to objectively quantify changes in brain states by comparing EEG signals before and after an intervention. The algorithm comprises two main steps: Spectral Exponent (SPEX) and generalized eigenvalue decomposition (GED). The SPEX component estimates the excitatory-inhibitory (E-I) balance in the brain by analyzing the power spectral density (PSD) of EEG signals, while the GED component identifies the optimal subset of EEG electrodes that reflects changes in the E-I balance before and after the intervention. 
The concept of SPEX was developed in silico in computational neuroscience (Gao et al., 2017) and validated using local field potential (LFP) of rat CA1 (Mizuseki et al., 2011), electrocorticogram (ECoG) of macaque monkeys (Yanagawa et al., 2013), and human EEG and MEG datasets (Donoghue et al., 2020). A voltage-clamp experiment determined that the time constant of excitatory glutaminergic AMPA receptors, τ_e, is 2.6-2.8 ms, while that of inhibitory GABAA receptors, τ_i, is 8.0-10.5 ms in cortical Layer VI, Va, and Vb of cat parietal cortex (Destexhe et al., 2001). Under the Fourier transformation, AMPA’s sharp spikes in the time domain correspond to flatter frequency spectra, whereas GABAA’s broad spikes correspond to steeper frequency spectra. The significance of the SPEX is that it connects high-level meso- and macro-scopic electrophysiological measures with low-level generative mechanisms of the E-I balance, which is the mechanism in which pathological hyperexcitability occurs in various neuropsychiatric conditions.

Materials and Methods
EEG preprocessinig
Use of the median statistics in the customized Welch’s method, known as robust spectral estimation(Melman & Victor, 2016), helps to reduce the burden and provides workaround to the notorious EEG cleaning issue. The arithmetic mean is sensitive to outliers, which can be problematic when dealing with signals contaminated by impulsive noise or other non-Gaussian disturbances. The median is a robust statistic that is less affected by outliers, making it a suitable alternative for averaging in such cases.

Proposed algorithm
After performing all necessary preprocessing steps on the EEG signals for each subject, the power spectral density (PSD) is computed using Welch’s method. The data are segmented into 1-second non-overlapping windows, and a Hamming window is applied to each segment. PSDs are then computed for each segment and averaged across all segments to obtain a final PSD estimate. This procedure is implemented using the ‘spectopo’ function in EEGLAB (Delorme & Makeig, 2004). When visualized on a logarithmic scale, the typical shape of the PSD follows a 1/f pattern with a noticeable bump. This pattern can be decomposed into two components: (1) the peaks on PSD, which is referred as the ‘periodic component’ in the original reference paper, which corresponds to neural oscillations in specific frequency bands; and (2) the spectral exponent (SPEX), which is referred to as a ‘aperiodic component’ which is also referred to as 1/f component. The SPEX reflects the balance between excitation and inhibition (E-I balance) in neural circuits. This E-I balance shapes the overall dynamics of neural activity and influences the structure of the power spectrum.
The SPEX is modeled using the following mathematical equation as a function of frequency f: 

L(f)=b-log⁡(k+f^x)
(Eq. 1)
where b is an offset that shifts the overall power spectrum, k is the knee parameter which adjusts the curvature of the shape, x is the exponent which controls the slope of the aperiodic decay, and log represents the base-10 logarithm operator. 
We employed FOOOF (Fitting Oscillations & One-Over-F) algorithm (Donoghue et al., 2020) to extract the exponent which allowed us to quantify the balance of excitation and inhibition in neural circuits. After computing the individual SPEX for each subject, we constitute a group-level matrix. Let X_pre and X_post be m by n matrices, where m is the number of electrodes and n is the number of subjects. We then find the eigenvalues and eigenvectors by solving the following generalized eigenvalue problem:

RWΛ=SW
(Eq. 2)
where

R=X_pre X_pre^T
(Eq. 3)

S=X_post X_post^T
(Eq. 4)

W represents the eigenvectors, and Λ represents the eigenvectors. Using GED, we obtain eigenvectors that maximize the variance ratio between the two covariance matrices: 

argmax_w  (w^T Sw)/(w^T Rw)
(Eq. 5)
This identifies the direction of maximal variance difference between the matrices, which serves as a spatial filter. By taking the first component, we can obtain a subset of electrodes that show the most significant changes in activity due to the intervention.

Results

Discussion

References
 
Delorme, A., & Makeig, S. (2004). EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics including independent component analysis. Journal of Neuroscience Methods, 134(1), 9–21. https://doi.org/10.1016/j.jneumeth.2003.10.009
Destexhe, A., Rudolph, M., Fellous, J. M., & Sejnowski, T. J. (2001). Fluctuating synaptic conductances recreate in vivo-like activity in neocortical neurons. Neuroscience, 107(1), 13–24. https://doi.org/10.1016/s0306-4522(01)00344-x
Donoghue, T., Haller, M., Peterson, E. J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A. H., Wallis, J. D., Knight, R. T., Shestyuk, A., & Voytek, B. (2020). Parameterizing neural power spectra into periodic and aperiodic components. Nature Neuroscience, 23(12), 1655–1665. https://doi.org/10.1038/s41593-020-00744-x
Gao, R., Peterson, E. J., & Voytek, B. (2017). Inferring synaptic excitation/inhibition balance from field potentials. Neuroimage, 158, 70–78. https://doi.org/10.1016/j.neuroimage.2017.06.078
Melman, T., & Victor, J. D. (2016). Robust power spectral estimation for EEG data. Journal of Neuroscience Methods, 268, 14–22. https://doi.org/10.1016/j.jneumeth.2016.04.015
Mizuseki, K., Diba, K., Pastalkova, E., & Buzsáki, G. (2011). Hippocampal CA1 pyramidal cells form functionally distinct sublayers. Nature Neuroscience, 14(9), 1174–1181. https://doi.org/10.1038/nn.2894
Yanagawa, T., Chao, Z. C., Hasegawa, N., & Fujii, N. (2013). Large-scale information flow in conscious and unconscious states: an ECoG study in monkeys. Plos One, 8(11), e80845. https://doi.org/10.1371/journal.pone.0080845 

