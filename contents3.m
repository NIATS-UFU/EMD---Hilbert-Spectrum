% Contents file
% by Adriano de Oliveira Andrade
% aoandrade@yahoo.com.br
% a.d.o.andrade@rdg.ac.uk
% Sep 2004
%
% Empirical Mode Decomposition Toolbox
%
% Remarks:
%         (1) In order to speed up the estimation of IMFs, DLLs were created from *.m files, 
%             therefore, any modification in *.m files will require a new compilation
%
% Function name                  Description 
%
% sig_to_imf                     It decomposes an input time-series x into its intrinsic mode functions
% imf_est                        Estimation of a single intrinsic mode function (IMF)
% meanEnv                        Estimation of the average envelope of a time-series
% zero_crossing.m                Estimation of the number of zero-crossings of an input sequence y
% sng.m                          Signal function
% extrema.m                      Estimation of local maxima (peaks) and minima (valleys) of time-series
%
% Demos
%
% EMD_demo1.m                    Illustration of the Empirical Mode Decomposition method
%
%filtLPD_demo1.m                 It shows how to estimate the frequency response of a low-pass differential filter  
%
% wavelet_demo1.m  
% wavelet_demo2.m                 This example illustrates the relationship of filters to 
%                                 wavelet shapes(required file: waveletDemo2.mat)
% wavelet_demo3.m 
% wavelet_demo4.m   
% wavelet_demo5.m   
%
% filtEMD_demo1.m
% filtEMD_demo2.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo3.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo4.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo5.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo6.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo7.m  (required functions: funcFiltEMD_demo2.m and funcFiltEMD2_demo2)
% filtEMD_demo8.m  (required file: filtEMD_demo8.mat)
%
% IMFsPhysicalMeaningDemo1
% IMFsPhysicalMeaningDemo2
% IMFsPhysicalMeaningDemo3
% IMFsPhysicalMeaningDemo4
% IMFsPhysicalMeaningDemo5
%
% End of Empirical Mode Decomposition Toolbox
%
%
% Hilbert Spectrum Toolbox
% Remarks:
%
% Function name                  Description 
%
% convScale.m                    Convert a scalar or a vector to a different scale
% DS.m                           Estimation of the degree of stationarity
% plotHS1.m                      Hilbert Spectrum generation 
% instantAtrib.m                 Instantaneous atributes of a time series
% MHSpec.m                       Marginal Hilbert Spectrum estimation
%
% End of Hilbert Spectrum Toolbox
%
%
% General functions 
% Remarks:
%
% Function name                  Description 
%
% histogram.m                    Estimation of the histogram of a time series
% mainPeaks.m                    This function estimates the principal peaks of an input time-series
% rmsEst.m                       root-mean square estimation of a time-series
% AR_LMS.m                       AR model estimation (based on the LMS algorithm - old version)
% arlms.m                        AR model estimation (based on the LMS algorithm) 
% spectro.m                      power spectrum estimation based on the AR model 
%
% End of General functions 
%
%
% EMG Signal Decomposition Toolbox
% Remarks:
%
% .............                  Signal thresholding and selection
% Function name                  Description 
%
% DigitalTrigger.m  
% threshold.m                     
% ss_peakDetector.m             detection of region of activities 
% ss_plotRAs.m                  it plots region of activities 
%
% ss_alignMUAPs.m             
% ss_get_alignedMUAPs.m          Extract raw and aligned raw MUAPs from signals based on interpolation 
% ss_get_alignedMUAPs2.m         Extract raw region of activities (RAs) and aligned raw RAs from signals        
% ss_get_alignedRA.m              
% ss_getRAs                      This function extracts region of activities from a given EMG signal. 
% ss_plotAlignedRAs              tool for visualization of region of activities
%
% ss_get_fixedSizeWindowFromRA.m This function extracts a fixed size window from region of activities 
% ss_wndCentredOnMainPeak.m      This function is similar to ss_get_fixedSizeWindowFromRA.m, but it extracts information
%                                from matrices and not directly from a time-series
%
%
% .............                  Signal Filtering
% Function name                  Description 
%
% ss_filtWavelet.m               Signal filtering based on the stationary wavelet and soft-threshold
% ss_filtWavelet2.m              Signal filtering based on reconstructed approximation and detail coefficients and soft-threshold
% ss_filtEMD.m                   Signal filtering based on the Empirical Mode Decomposition method
% ss_filtDiff.m                  Differential filter
% SigDenoise.m                   Signal filtering based on the Empirical Mode Decomposition method (old)
% ss_backgroundActRemoval.m    
%
% .............                  Feature extraction
% Function name                  Description 
%
% PCAanalysis.m                  GUI for PCA estimation
% ss_timeDomainFeatures.m        features based on the time-domain properties of MUAPs
% ss_FrequencyDomainFeatures.m  
% ss_getContWavFromRefVec.m       
% ss_threshold_crossing.m         
% ss_estRiseTime.m               
%
% .............                  Clustering
% Function name                  Description 
% ss_clusterEMG                  This function cluster EMG signals based on GTM (automatic method)
% GTMscript.m 
% GTMscript1.m          
% GTMscript2.m      
% ss_TrainSOM.m    
% ss_trainGTM.m                   this function fits a GTM model to a given data set 
% ss_initGTM_neuralGas.m             
% ss_DefineBoundSOM.m      
% ss_labelDataSetGTM.m            it assigns patterns to latent points/nodes based on probabilities   
% ss_assignPatternToGTMCluster.m
% ss_plotLSgrid.m                 
% ss_contiguityGTM.m            
% ss_ExtracRefVecFromGTM.m     
% ss_buildGTMclustering.m      
% ss_showLatGrid.m                
% ss_SOMclusterSim.m         
% ss_EstimateMeanPattern.m   
% ss_getMUAPtemplates.m           
% ss_EstiamteMeanPattern.m    
% ss_classifyRefVecs.m          
%              
% ss_getMUAP.m      ??????????              
% ss_gsamp.m        ??????????     
%
%       
% ss_getRA_ARLMS.m            ??????????    
% ss_ContWaveletFeatures.m    ?????????
% ss_getRA_IPI.m              ?????????     
% ss_get_ARcoefMUAPs.m        ????????     
% ss_get_RA_splines.m         ???????    
% ss_ExtractInterpFeatures.m  ???????  
%         
% ss_normal.m 
% ss_ICA.m  
% ss_ppNEM.m     
% ss_featExtraction_1.m    
% ss_sincFilter.m       
% ss_mixPat.m         
% ss_powerSpectrumMUAPs.m    
%       
% ss_detectOverlap.m ?????????????????????           
%
% Demos                           Description 
%
% ss_demo1.m                      it shows how GTM clusters data from 10 different set of MUAPs
% ss_demo2.m                      selection of region of activities, threshold estimation and problems regarding the peak detector 
% ss_demo3.m                      
% ss_demo4.m  
% ss_demo5.m                      it shows how GTM deals with experimental EMG signals
% ss_demo6.m  
% ss_demo7.m  
% ss_demo8.m  
% ss_demo9.m  
% ss_demo10.m  
% ss_demo11.m                     Analysis of some experimental data from Japan
%
% End of EMG Signal Decomposition Toolbox
%
%
%
% EMG simulation toolbox
% Remarks: (1) Two different approaches are provided by this toolbox. The first one is the generation of
%          the EMG signal based on the single fibre action potential of Rosenfalck and on the 
%          very well known structural model for EMG generation. The second approach is the simulation
%          of EMG signal based on experimental data. For this, single MUAPs and windows of noise
%          were selected from different sets of EMG signals. These data are used to generate EMG signals
%          based again on the very well known structural EMG model
%          (2) The best way of understanding how the functions bellow work is to study how they
%          are called in the main graphic user interface (guidsp.m)
%
% EMG model using experimental data
% Function name                  Description 
%
% exp_GenerateEMG.m              Generation of EMG signals based on experimental data
% exp_GenerateMUAPT.m            Generation of MUAPTs based on experimental data
% getPatternFromMUAPs.m          Generation of Patterns of MUAPs based on the Point Distribution Model 
%
% EMG model using solely synthetic signals
% Function name                  Description 
%
%..............                  Fibre action potential 
%
% EAP.m                          ???
% extPot.m                       Simualtion of the single fibre action potential  
% GenerateFAP.m                  Generation of single fibre action potential (with GUI interface; it uses extPot.m)
% GUI_AP.m                       GUI for single fibre action potential simulation
% siap.m                         GUI for single fibre action potential simulation (old)
%
%
%..............                  Muscle cross-section and electrode positioning
%
% EstFibreToElecDist.m           Estimation of distance from electrodes to fibres 
% GenerateCircle.m               Generation of a circle at arbitrary position with a given diameter
% GenerateMuscle.m               Generation of muscle cross-section       
% SetElectrodePos.m              Electrode positioning
% PointIsInTheCirle.m            Verifies whether a point is in the circle or not
%
%..............                  Motor unit action potential
%
% muap.m                         Motor unit action potential simulation
% SimulateMU.m                   Simulation of muscle cross-section, electrode positioning and MUAP templates
% GenerateMUAPpatterns.m         This function generates a set of MUAPs based on a single MUAP template
%                                (It uses the Power spectrum of the signal for this purpose)
% ScaleMUAP.m                    Generation of Patterns of MUAPs based on MUAP template (Using signal power spectrum)
% VisualizeMUAPtemplates.m       Visualization of motor unit action potential templates
% GUI_MUAP.m                     GUI for motor unit action potential generation
%
%
%..............                  Motor unit action potential train
%
% GenerateMUAPT.m                Generation of motor unit action potentials (GUI, it uses muapt.m)
% muapt.m                        Generation of motor unit action potentials
%
%..............                  Motor unit action potential train
%
% VisualizeMUAPT.m               Visualization of MUAPTs
% GUI_MUAPT.m                    GUI for motor unit action potential train
% VisualizeMUAPTfiringTime.m     Tool for visualizing MUAPT firing time
% MUAPT_IPIhist.m                Tool for visualizing MUAPT inter-pulse interval histogram
%
%..............                  Composite EMG signal
%
% SimEMG.m                       Simulation of the composite EMG signal
% SimulateCompEMGsignal.m        Simulation of the composite EMG signal
% GUI_EMG.m                      GUI for composite EMG signal simulation
% guisemg.m                      GUI for composite EMG signal simulation (old)
%
%..............                  Utilities 
%
% OverlapBuilder.m               Tool for building MUAP overlaps
% MUAPoverlapBuilder.m           Tool for building MUAP overlaps
% MUAPOverlapsBuilder.m          Tool for building MUAP overlaps
% gauss.m                        Estimation of a normalized Gauss distribution. 
% sim_labelRA.m                  This function labels region of activities (RAs) of
%                                synthetic EMG signals
% sim_labelRA2.m                 This function labels sub-region of activities of
%                                synthetic EMG signals
%
% Demos
%
% simEMGdemo1.m    
% simEMGdemo2.m    
%        
% End of EMG simulation toolbox
%
%
% GMM (Gaussian Mixture Model) toolbox
% Function name                  Description 
%
% gmm_EM.m                       This function implements the EM algorithm for GMM models
% gmm_dist.m                     This function estimates the squared Euclidian distance between 2 matrices 
% gmm_prob.m                     This function computes the data probability p(x)
% gmm_activ.m                    This function computes the component likelihoods (or ativations)
% gmm_kmeans.m                   This function implements the k-means algorithm
% gmm_create.m                   It creates a Gaussian mixture model
% gmm_post.m                     This function computes the component posterior probabilities
%
% Demos                       
%
% gmm_demo1.m                    
% gmm_demo2.m                    
% gmm_demo3.m  
% gmm_demo4.m    
%
% End of GMM (Gaussian Mixture Model) toolbox
%
%
% Digital Signal Processing Toolbox (Graphic User Interface)
% Remarks: (1) The graphic user interface for this toolbox is implemented in the files guidsp.m
%              and guidsp.fig 
%          (2) Only the selected signal is processed. Multiple selection is not allowed
%
% Menu Basic Tools               Description 
%
% Extrema                        Estimation of the extrema of the current signal by using the function 'extrema.m'
% Histogram                      Histogram of the current signal ('histogram.m' is used)
% Envelope                       Estimation of signal envelopes. The function 'meanEnv.m' is employed
% Burst                          Burst detection (detectBurst.m)
% Principal Peaks                Estimation of the main peaks of a time-series (mainPeaks.m)
% Instantaneous atributes        Estimation of the instantaneous atributes of a time-series (instantAtrib.m) 
% RMS                            Root-mean square of a time-series (rmsEst.m)
% Power spectrum                 Fourier power spectrum of a signal (the function PSD - Matlab - is used)
% Spectrogram                    Estimation of the spectrogram of a signal (the function specgram - Matlab - is used)
%
%
% Menu Hilbert Analysis          Description
%
% Empirical mode decomposition   Estimation of intrinsic mode functions (sig_to_imf.m)
% Hilbert spectrum: 
%             2-D plot           Hilbert spectrum estimation (plotHS1.m, type = 0)
%             2-D contour plot   Hilbert spectrum estimation (plotHS1.m, type = 1)
%             3-D contour plot   Hilbert spectrum estimation (plotHS1.m, type = 2)
% Marginal Hilbert Spectrum      Estimation of the MHS (MHSpec.m)
% Degree of stationarity         Estimation of the degree of stationarity (DS.m)
% IMF analysis                   Power spectrum of different IMFs
% Completeness analysis          Reconstruction of a signal based on the summation of IMFs
%
%
% Menu EMG simulation            Description
%
% External fibre AP              Simulation of muscle fibre action potential
% Motor unit AP:
%         MUAP simulation        Simulation of motor unit action 
%         MUAP visualization     Visualization of generated motor unit action potentials
%         Pattern generation     Generation of Patterns of MUAPs based on MUAP template (ScaleMUAP.m)
% Motor unit AP train:
%         MUAPT simulation       simulation of motor unit action potential train
%         MUAPT visualization    visualization of simulated motor unit action potential train
%         Visualize firing time  visualization of the firing time of MUAPTs
%         IPI histogram          inter-pulse interval histogram
% Composite EMG                  simulation of the composite EMG signal
% Muscle cross-section           visualization of the muscle cross-section 
% MUAP overlap builder           tool for simulating MUAP overlaps
% Data-driven model:
%         Generate MUAPT         MUAPT simulation based on experimental data
%         Generate sEMG          simulation of synthetic EMG signals based on experimental data  
%
%
% End of Digital Signal Processing Toolbox (Graphic User Interface)
%            
% Butfilter.m                    Low pass Butterworth filter 
% changeVec.m                    NOTHING
% detectBurst.m                  NOTHING
%
% envelope.m                     Estimation of the upper and lower envelopes of time-series
% EnvMAF.m                       ?????
% EstIMFsBatch.m                 ?????
% filter.m                       ?????
% FresponseWLPD.m                ?????
% IPIplot.m                      ???
% LoadDataFile.m                 ???
% LoadFile.m                     It loads a text file
% LowPassBut.m                   ???
% MAF.m                          Moving average filter
% meanEnvBackUP.m                NOTHING
% meanEnvelope.m                 ????
% meanEnvelope1.m                ????
% movavgfilter.m                 ????
% plotHS.m                       It generates the Hilbert Spectrum of a signal
% readDSCfile.m                  ????
% RMS.m                          ????
% SearchList.m                   ????
% signum.m                       ????
% str2cell.m                     ???
% threshold_crossing.m           ???
% thresholdOld.m                 ???
% WLPD.m                         weighted low pass differential filter