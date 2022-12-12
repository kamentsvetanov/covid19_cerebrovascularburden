function [S,T,errors] = wtcvr_04_postproc_func(S,T,dataname,wave,signalType,doRun)
% Function to post-process EPI data and generate RSFA maps


%% 3. Post-Process EPI data
P_name = {'Name'        ...  1. Name for processing combination
         'detrend'      ...  2. Add detrend
         'regressCovs'  ...  3. Regress Covariates
         'regressGlobal'...  4. Regress Global signal
         'regressWM'    ...  5. Regress WM signal
         'regressCSF'   ...  6. Regress CSF signal
         'CompCor'      ...  7. Regress CompCorr Components (Behzadi et al 2007)
         'derivatives'  ...  8. Add derivatives
         'squaredTerms' ...  9. Add squared terms
         'motionOpt',   ... 10. Motion option 0 - No motion;  1- Power et al'13; 2 - Satthertwaite et al '12 (24 parameters in GLM); 3 - 6 RPs
         'filterMethod' ... 11. Add filtering 0 - No filtering;  1 - Do filtering part of regression (Hallquist et al 2013) | 2 - Butterworth filter, sequential regression
         'filterType'   ... 12. 1 - Band-pass; 2 - High-pass; 3 - Low-pass 
         'PreWhiten'    ... 13. PreWhiten
         'zscoreY'      ... 14. Zscore input EPI timeseries before processing
         'zscoreYr'     ... 15. Zscore output EPI timeseries after processing
         'zsoreCovs'    ... 16. Zscore covariates/all regressors
         'wantZval'     ... 17. Zval connectivity
         'trimVolumes'  ... 18. Trim N initial volumesn
         'addConstant'  ... 19. Add constatnt value to all voxels to hack 1st level SPM
         'saveSD'       ... 20. Compute temporal SD of timeseries (e.g. for QC) or RSFA maps
         'saveTS'       ... 21. Save processed timeseries (or rather not if interested only in RSFA maps)
         'saveFC'       ... 22. Compute Functional connectivity between ROIs
         'saveCVR'      ... 23. Compute rCVR as per Liu et al 2017 Neuroimage
         'HPC'          ... 24. High-pass filter
         'LPC'          ... 25. Low-pass filter
         };


%    1.              2     3   4   5   6    7    8    9    10       11       12   13   14     15       16         17      18       19      20      21    22       23   24   25                        
%   Name           Detr regCov GS WM  CSF  Cc  Deriv ^2  movOpt FiltMeth  FiltTyp PW zscoreY zscoreYr zsoreCovs wantZval trimVol addConst saveSD saveTS saveFC saveCVR HPC  LPC   
P ={
    'wcm2f11_10_128' 1     1    0  1   1   1     1    1    2         1       1     0    0      0        1          0       5        0       0      1      0       1   .0078 .1 ;   ... Processing for Alistair
 };

P = cell2table(P,'VariableNames',P_name);

for pp = 1:size(P,1)%numel(postproc_name)
    S.process.postproc(pp).name               = P.Name{pp};
    S.process.postproc(pp).do                 = table2struct(P(pp,:));
    S.process.postproc(pp).LPC                = 1/P.LPC(pp);
    S.process.postproc(pp).HPC                = 1/P.HPC(pp);
    S.process.postproc(pp).TR                 = S.data.restingW1.TR;
    S.process.postproc(pp).smooth             = [6 8 12];
end

numProcess  = size(P,1);
Ns          = size(T,1);
fpostproc   = cell(size(T,1),2);

% ---------------------------
% Predefining some variables
% ---------------------------

nameData    = [dataname 'W' num2str(wave)];  
rootDir     = sprintf('%s/restingW%d/data_postproc_%s/',S.paths.postproc,wave,lower(dataname));%fullfile(S.paths.mriana,nameMeasure,['data_rsfa' num2str(iWave) '_' lower(nameData)]);
f_nameEPI   = ['f_resting' dataname 'W' num2str(wave)];
f_nameSignal = ['f_node_resting' dataname 'W' num2str(wave) '_hcp360'];
f_nameRP    = ['f_rpRestW' num2str(wave)];
f_mask      = S.process.f_brainMask;
f_maskGM    = S.process.f_gmMask;

errors = [];

for xx = 1:numProcess
    nameProc = S.process.postproc(xx).name;
    f_nameOut   = ['f_resting' dataname 'W' num2str(wave) '_' nameProc];
    % Loop through subjects
    parfor iSub = 1:Ns 
    % try 
        f_signal = T.(f_nameSignal){iSub};
        f_epi    = T.(f_nameEPI){iSub};
        if exist(f_signal,'file')
            [p fn e]  = fileparts(f_epi); 
            f_out = fullfile(rootDir,T.SubID{iSub},[nameProc '_p' fn '.nii']);     % Filtered and covaried out
            [p fn e]  = fileparts(f_out); 

            if ~exist(fileparts(f_out),'dir') | doRun
                mkdir(p);
                % Load Node and Tissue-based signals
                Signal   = load(f_signal);
                ySignal  = f_epi;%[Signal.ROI.(signalType)]; 
                yWM      = Signal.Tissue.wm.meanSignal;
                yCSF     = Signal.Tissue.csf.meanSignal;
                yCompcor = [Signal.Tissue.wm.pcaSignal Signal.Tissue.csf.pcaSignal];

    %             for xx = 1:numProcess
                    X        = [];
                    X        = S.process.postproc(xx);
                    X.dataset= nameData;
                    X.SubID  = T.SubID{iSub};
                    X.Y      = f_epi;% Input of EPI data
                    X.f_out  = f_out;
                    X.f_rp   = T.(f_nameRP){iSub};
                    X.f_in   = T.(f_nameEPI){iSub};
                    X.f_mask = f_mask;
                    X.f_maskGM = f_mask;
                    X.compcor= yCompcor;
                    X.csf    = yCSF;
                    X.wm     = yWM;
                    kat_fmri_postprocessing_GLM(X);
            end
            f_out_all{iSub} = f_out;
        end           
    end

    try
        T.(f_nameOut) = f_out_all';
    end
end
