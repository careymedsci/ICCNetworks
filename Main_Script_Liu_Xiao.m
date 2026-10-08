 
% ============================================================================
% Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
%
% This code was jointly written by:
% Thomas Broggini (Frankfurt University, Germany)
% Liu Xiao (Xiangyang First Hospital, China)
%
% The way how the mathematical models and methodologies implemented in this code have been
% individually customized and are not intended for generic use. 
%
% For permission requests or inquiries, please contact:
% Thomas Broggini  : broggini@med.uni-frankfurt.de
% Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
%
% Unauthorized use will be considered a violation of intellectual property rights.
% ============================================================================



% 作者：Thomas Broggini 和 刘晓； 法兰克福大学医院；湖北医药学院襄阳市第一人民医院
% By Thomas Broggini (Supervisor) and Xiao Liu (MD Student)
% Department of Neurosurgery
% Neuroscience Centre
% University Hospital Frankfurt
% Goethe University Frankfurt, Germany
% Frankfurt Cancer Institute, Germany

% Xiangyang No.1 people's Hospital, China
% Also an python version coded 
% xiao.liu@stud.uni-frankfurt.de
% 2024.10.09

% side_project_folder/
%     ├── main_analysis.m
%     └── function_module/ 
%         └──folders
%            ├── getAnalysisParameters.m
%            ├── loadROIData.m
%            ├── processPeaks.m
%             ...
%

% As I am not a professional programmer or data scientist, the code 
% formatting may be a bit messy，we appreciate your understanding

clc ; clear
%% ######################  Set working folder  ############################


% Navigate to the folder where the main function is located first
% Note! every moulde need these path settings!
mainPath = pwd;
functionPath = fullfile(mainPath, 'function_module');
if ~exist(functionPath, 'dir')
    error('Cannot find function_module folder at: %s', functionPath);
end
addpath(genpath(functionPath));  


%% #######################  Preparation works #############################


% Option 1: merge tiff stacks and extract red channel
merge_red_channel_stacks();
% Option 2: load merged stack
im_data = load_tiff_stack();
% Option 3: offset
%im_data = apply_offset_correction(im_data);
% Option 4: dF/F 
calculate_dfof_and_save(im_data, pwd);  
clc, clear all；


%% ####################  Figure 1 Analysis  ############################### 


% This code works for Figure 1
% Load tiff_stack first please
result_one_anaylsis(im_data);


%% ####################  Figure 2 Plotting  ############################### 


% only for Figure 2A 
Fig_2A_summary_cells_postive_ca_events_matalab()
% only for Figure 2B
Fig_2B_analyze_temporal_modes()


%% ################## Figure 3: Global Charactors analysis  ###############   


% Characteristics of intrinsically rhythmic Ca²⁺ activity globally
calcium_fft_analysis();

% ploting
coherence_power_analysis();

% Estimate significant coherence threshold (coherence error evaluation)； 
% Any peak of the coherence below this value cannot be picked up (Fig. 3B)
% 500 sims is the lowest requierment
coh_thresh = estimate_coherence_threshold(Vn, scores, tapers, Fs, nfft, ...
                                 f, findx, toplot.f_vector, 500, 0.05, 10);


%% ######################  Network Analysis  #############################


% ========================== step 1 do segment by Fiji (do not use imageJ)


%  if you cannot run it soomthly, go to Fiji run the macro independently
%  original Macro sciript named "new_fiji_processing.ijm"
[file_name, file_path] = uigetfile('*.tif', 'Select a video file');
full_file_path = fullfile(file_path, file_name);
disp(['Selected file: ', full_file_path]);
fiji_path = '"C:\Program Files\Fiji.app\fiji-windows-x64.exe"';
macro_file = 'fine_fiji_processing.ijm';
full_macro_path = fullfile(functionPath, filesep, macro_file);
full_macro_path = ['"' full_macro_path '"'];
ij_args = ['inputdir=' file_path ',filename=' file_name];
command = [fiji_path ' -macro ' full_macro_path ' "' ij_args '"'];
system(command);


% ========================          step 2  set parameters       =========


clc
clear
% Reminder: Model 1, 2, 3 are only controls used to evaluate the presence of 
% true calcium signaling in experimental data correlation analysis (step 4). 
% They should not be used for global analysis, as any recombination of the 
% real data will inevitably result in some false positive correlations 
% between some ROIs, leading to potential false positives.
% Only you have to caculate the long distance ROIs, model 4 can be chose to
% do the evaluation (like you have to set above 150um to fit the powerlaw and so on)
% Please check the mauscript to understand the parameters and preset them
% suitable.
params = getAnalysisParameters();

% loading the rois
% read file named: "xxxx_mean_values.csv or xxxx_mean_values.txt" which 
% was gerated by Fiji
params = loadROIData(params);


% =========================        step 3  process the peaks     ==========


[params, analysis_summary] = processPeaks(params);

% Optional, choose the plot you want===
% load('all_ver_0.mat');
% % plot traces simple version
% plotSimpleTraces(params);
% % plot traces fine version
% plotEnhancedTraces(params);
% % visulize the selected cells
% plotChosenCellsonMaxZimage(params);


% =========================        step 4  Correlations caculation   ======

dcutoff_final = compute_local_dcutoff(rparams, params);
% This is only for distance cutoff evaluation; But just a reference.

% load('all_ver_0.mat');
[params, analysis_summary] = processCorrelations(params, analysis_summary); 

% small-world parameters caculation （C/L caculation）
[params, analysis_summary] = network_theory_analysis(params, analysis_summary);
      
 % KS test needed 
 % Please run all the modes except 4 first!!!
 % please save the picure manually!
 % For the Figure. 4C
 % KS_CDFs_test(params);


% =========================        step 5  Visualization      ============


% **********    Cutoff evaluation and  Connectivity Analysis   ************

% CRITICAL! BEFORE VISUALIZING THE CORRELATION, PLEASE CHECK THE CORRELATION 
% THORESHOLD (CUTOFF) PRESET. IF IT IS ABOVE THE HIGHEST CUTOFF CALCULATED 
% BY MODE 1-3, THEN CONTINUE WITHOUT ANY CHANGE; OTHERWISE, CHANGE To
% THE CUTOFF BY SELECTING THE HIGHEST VALUE.

% Way to check the cutoff: check these figures
% 1, * _circular_shift_corr_vs_distance1.fig;
% 2, * _scrambled_corr_vs_distance2.fig;
% 3, * _linear_shift_corr_vs_distance3.fig
% If these marked value both less than what you have preset by parameters
% "params.threshcorr", then just countinue without any change or change to
% the highest value you get from the 3 figures.
[params, analysis_summary] = analyzeNetworkConnectivity(params, analysis_summary);


% *******  Visualization the correlations and network components **********


% change the last number according to the
% meaning of the mode. '0' is the empirical data!

% load('all_verVisulization_0.mat'); % if you lost you data from the workspace 
testback = 'no';
params.testback = testback;
[params, analysis_summary] = visualizeCorrelations(params, analysis_summary);
% % Visualization == hub cells connctions （optional, just check the topology）
% [params, analysis_summary] = visualizeCorrelations_hub_connections(params, analysis_summary);
% % Visualization ≥ hub cells connctions （optional, just check the topology）
% [params, analysis_summary] = visualizeCorrelations_hub_above(params, analysis_summary);


% first caculation, please just do it. This is for
% mix-mFFT-DPSSVD-corss-correlation analysis below
testback = 'no'; 
params.testback = testback;


% Find the trigger cells（assumed settings and definitions, for reference only）
[params, analysis_summary] = findTriggerCells(params, analysis_summary);

% Find the hub cells
[params, analysis_summary] = findHubCells(params, analysis_summary);


% =========================      step 7  Periodicity theory and visualization
testback = 'no';
params.testback = testback;
[params, analysis_summary] = findPeriodicCells(params, analysis_summary);
[params, analysis_summary] = analyze_peakwidths_of_p_and_active(params, analysis_summary);


% =========================        step 8  plot 3D traces    =============


% Plese run all the code before this sentence
% Please note that here all the figures will not save itself, please 
% save your necessary figures
load('all_for_3D_0.mat');

% if there is error showing cannot find 'all_ver_0.mat', please check the
% code and correct the name to 'all_ver.mat'
[params, analysis_summary] = visualize3DCellTraces(params, analysis_summary);


% =========================        step 9  Classifying cells   ==========
% based on the number of neighbors
[params, analysis_summary] = classifyCellsByNeighbors(params, analysis_summary);


% =========================        step 10  Identifying small-world   ====
% characteristics by mathematical principles
%  Please note that here all the figures will not save itself, please 
%  save your necessary figures


% Fig.14F (only sigle data for evaluation, different with plot_k_vs_sigma()) 
% Fig.15A 
% set preset parameter xmin = 1 please;  
load("all_for_analysis_0.mat")
[params, analysis_summary] = plotPowerLawDistribution(params, analysis_summary);

% degree distribution (use C matrix), they should have same relusts
% Fig.15A (only sigle data for evaluation）
plotDegreeDistribution(params); 

% Fig. 15B
% This is a code for scale free and small world test to 
% all biological replicates and technical replicates
% small world visualization
% Before run codes below in this paragraph, please prepare all the data
S = plotSmallWorldMetrics(); % different with small_word_test();

% Fig. 14F
% Before run codes below in this paragraph, please prepare all the data
% here we analysis by file '* _probability_distribution.txt'
% ofcourse you can also use the I already precaculated directly
% different with plot_k_vs_sigma_drug() below;
plot_k_vs_sigma(); 


% ===================  step 11  network weankness analysis  =============
% Fig. 16E and F
[params, analysis_summary] = weaknessNetwork(params, analysis_summary); 


%% ###########################  layer analysis  ###########################

% Note: Before run this, please prepare all the data
% figure 16B， read Excel data file name: "summary_on_all_math_fit.xlsx"
periodicity_ratio_Hub_vs_non_hub(); 

% Fig.16H
% different from "plot_freq_Min_Max_group_comparison()"
[params, analysis_summary] = plotFrequencyHistogram(params, analysis_summary);

% for data colloction;
% file name: summary_on_all_math_fit.xlsx
exportMathFitToExcel();

% Fig.16G
% read Excel data file name: "summary_on_all_math_fit.xlsx"
% this file will be colloct later, do this code after all the data
% processed
directionaryscatterWithErrorBars();


%%  ############  mix mFFT with cross-correlation analysis  ############### 


% step 1 get cell mfft pifctures


% figure 17
calcium_fft_analysis_combine();
% load all results from correlation matrix
load('all_verVisulization_0.mat');

testback = 'yes';
params.testback = testback;
[params, analysis_summary] = visualizeCorrelations(params, analysis_summary);
% % Visualize == hub cells connctions （optional, just check the topology）
% [params, analysis_summary] = visualizeCorrelations_hub_connections(params, analysis_summary);
% % Visualize ≥ hub cells connctions （optional, just check the topology）
% [params, analysis_summary] = visualizeCorrelations_hub_above(params, analysis_summary);

% do not skip analysis trigger cells, even you dont care about it
[params, analysis_summary] = findTriggerCells(params, analysis_summary);
[params, analysis_summary] = findHubCells(params, analysis_summary);
[params, ~] = findPeriodicCells(params, analysis_summary);


% step 2  run all the code above in this part beforehand


 analyze_frequency_distribution_mFFT();
 % Fig. 17C and D
 analyze_hub_power_distribution(toplot, params, x_coords, y_coords); 
 % Fig. 17H and I
 analyze_periodic_power_distribution(toplot, params, x_coords, y_coords);
 % Fig. 17E and F


 % step 3 summary


 % Figure 18A and B
 % Check the direction of lines on correlation first
 % Please note that here all the figures will not save itself, please 
 % save your necessary figures
 load('all_verVisulization_0.mat');
 testback = 'no';
 params.testback = testback;
 visualizeCorrelationsWithArrows(params);

 % plot the cells you want
 % Fig. 18C
 plot_cell_traces(all_traces, 2);

 % plot the network layer
 % read Excel data file name: "summary_on_all_math_fit.xlsx"
 layer_analysis();


%%  ################  treatments analysis and others  ##################### 


% Fig. 19---  reagents treatment analysis 
% 300 frames comparation with baseline
params = getAnalysisParameters();
params = loadROIData(params);
[params] = processPeaks(params);

% 3D plots for control group; eg. Fig. 19B 
plot3DTraces_r(params);
% 3D plots for treatment group； eg. Fig. 19C
plot3DTraces_r_t(params);
% for positive reagents ； eg. Fig. 19C
plot3D_positive_reagents(params, savepath);
% pick randomly traces for 3D plots to highlight the overall shape
plot_csv_waveforms();


% 200 frames for comparsion analysis, note that some title are changed 
% and peaks of ROIs defined as at leaset 3 peaks
% please read the manuscript first
params = getAnalysisParameters();
params = load200ROIData(params);
[params, analysis_summary] = processPeaks(params);
[params, analysis_summary] = processCorrelations(params, analysis_summary); 
[params, analysis_summary] = network_theory_analysis(params, analysis_summary);   
[params, analysis_summary] = analyzeNetworkConnectivity(params, analysis_summary);
testback = 'no';
params.testback = testback;
[params, analysis_summary] = visualizeCorrelations(params, analysis_summary);
[params, analysis_summary] = findTriggerCells(params, analysis_summary);
[params, analysis_summary] = findHubCells(params, analysis_summary);
[params, analysis_summary] = findPeriodicCells(params, analysis_summary);
[params, analysis_summary] = classifyCellsByNeighbors(params, analysis_summary);
[params, analysis_summary] = plotFrequencyHistogram(params, analysis_summary);
[params, analysis_summary] = analyze_peakwidths_of_p_and_active(params, analysis_summary);
[params, analysis_summary] = plotPowerLawDistribution(params, analysis_summary);
save('_all_result', 'params', 'analysis_summary');


%%%%%%%%%%%%%%%
% for data collection;
% file name: pharmacological investigation.xlsx
exportAnalysisToExcel();

%%%%%%%%%%%%%%%
% for data colloction;
% file name: summary_on_all_math_fit.xlsx
exportMathFitToExcel();

% Comparsion results;
% Statstical results name: treatment_comparison_stats.mat/xlsx
% eg. Fig19D-L and N 
results_treatemnt_stats_ = drugTreatmentsub();

% for B16-F10.n3; file name:co-culture_in_vivo_characteristics.xlsx
peaks_and_coordinating(); 

% file name: pharmacological investigation.xlsx
% eg. Fig19M 
% Statstical results name: freq_results_stats.mat or 
% frequency statistical_results.xlsx
results_freq_stats = plot_freq_Min_Max_group_comparison();

% file name: pharmacological investigation.xlsx
% eg. Fig19O and S
plot_k_vs_sigma_drug();

% file name: pharmacological investigation.xlsx
% eg. Fig19R 
small_word_test();

% file name: pharmacological investigation.xlsx
% eg. Fig19P and Q
periodicity_ratio_Hub_vs_non_hub_drug();

%(mind the row and col)
% co-cultrure and in vivo analysis
% eg. Fig31
% Statstical results name: KW_results.mat
results1 = analyze_co_culture();
% eg. Fig33; 40
% Statstical results name: KW_results.mat
results2 = analyze_animal_part();
