function [S,T] = fmricv_02_process_resting(S,T)


%% Smooth Normalized data by 6mm (required for ICA-AROMA)


smooth_filter = 6; % Smoothing filter in mm (6mm as per Pruim et al 2015 paper)
S.process.aroma.smooth_filter = 6;

Waves = [1];

for iWave = Waves
    f_smooth = cell(size(T,1),1);
    nameMeasure     = ['restingW' num2str(iWave)];
    nameData        = ['f_restingNormW' num2str(iWave)];
    
   parfor iSub = 1:size(T,1)
      
        f_in = T.(nameData){iSub};
        [p,f,e] = fileparts(f_in);
        if exist(f_in,'file')
              
            fn = fullfile(S.paths.postproc,nameMeasure,'data_smooth',T.SubID{iSub},sprintf('s%d%s%s',smooth_filter,f,e));
            f_smooth{iSub} = fn;
            
            if ~exist(fn,'file')
            
                mkdir(fileparts(fn));


                % -------------------------------------------------------------
                % Normalise the entire data set (the average across all voxels 
                % and all time points) is being scaled so that that value is 
                % set to 10000 (FSL's 4D global normalisation)
                % -------------------------------------------------------------
                V = spm_vol(f_in);
                Y = spm_read_vols(V);
                normFactor = 10000/mean(Y(:));
                Ynorm = Y.*normFactor;

                % Save 4D image
                ind  = cat(1,V.n);
                N    = cat(1,V.private);

                mx   = -Inf;
                mn   = Inf;
                descr= [];
                for i=1:numel(V),
                    str1     = spm_str_manip(V(i).fname,'rt');
                    descrip  = str1(16:end);
                    descr    = [descr ' ' descrip];
                end

                mx         = max(Ynorm(:));
                mn         = min(Ynorm(:));
                sf         = max(mx,-mn)/32767;
                ni         = nifti;
                ni.dat     = file_array(sprintf('%s',fn),[V(1).dim numel(V)],'INT16-BE',0,sf,0);
                ni.mat     = N(1).mat;
                ni.mat0    = N(1).mat;

                ni.descrip = descr;%'4D image';%
                create(ni);
                for i=1:size(ni.dat,4),
                    ni.dat(:,:,:,i) = Ynorm(:,:,:,i);% N(i).dat(:,:,:,ind(i,1),ind(i,2));
                    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
                end

                % ----------------------
                % Smooth image
                % ----------------------
                spm_smooth(fn,fn,repmat(smooth_filter,3,1));
            end
        end
    end
    T.(['f_restingSmoothW' num2str(iWave)]) = f_smooth;
end


%% 3. Run ICA-AROMA 
% Preprocessing steps required for date before applying AROMA as per Prium et al 2015 Neuroimage is:
% - motion correction
% - normalization
% - smoothing 6mm


% --------------
% Aroma settings
% --------------
S.process.aroma.fnames.aromapy  = '/imaging/camcan/sandbox/kt03/other/toolboxes/ICA-AROMA-master/ICA_AROMA.py';
S.process.aroma.brain_mask      = S.process.f_brainMask;%'/imaging/camcan/sandbox/kt03/fmri/masks/BrainMask_61x73x61.nii';% /imaging/camcan/sandbox/ds02/templates/masks/brain/mask_ICV_61x73x61.nii'; % Mask with a subjects full brain coverage. Ideally would be subject specfic BET output
S.process.aroma.num_ses         = 1;


for iWave = Waves
    
    f_aroma         = cell(size(T,1),1);
    dir_aroma_add   = cell(size(T,1),1);
    nameMeasure     = ['restingW' num2str(iWave)];
%     if strcmp('release002',S.project.release)
        nameData        = ['f_restingSmoothW' num2str(iWave)]; %'f_restingSmoothW' - release002;  'f_restingNormW' - release003
%     else
%         nameData        = ['f_restingNormW' num2str(iWave)]; %'f_restingSmoothW' - release002;  'f_restingNormW' - release003
%     end
    nameRP          = ['f_rpRestW' num2str(iWave)];
    nameAroma       = ['f_restingAromaW' num2str(iWave)];
    rootDir         = fullfile(S.paths.postproc,nameMeasure,'data_aroma');
    
    parfor iSub = 1:size(T,1)
            
        f_in = T.(nameData){iSub};
        
        if exist(f_in,'file')
            [~,fn]              = fileparts(f_in);
            f_aroma{iSub}       = fullfile(rootDir,T.SubID{iSub},['d' fn '.nii']);
            
            % If file doesnt exist run aroma
            if ~exist(f_aroma{iSub},'file')
                dir_aroma_add{iSub} = fullfile(fileparts(f_aroma{iSub}),'additional/');
                mkdir(dir_aroma_add{iSub});
                
                X.f_in         = f_in;
                X.f_out        = f_aroma{iSub};
                X.dir_aout     = dir_aroma_add{iSub};
                X.f_rp         = T.(nameRP){iSub};
                X.aroma_path   = S.process.aroma.fnames.aromapy;
                X.f_brain_mask = S.process.aroma.brain_mask;
                X.TR           = S.data.(nameMeasure).TR;
                ds_run_aroma(X);
            end
            X = [];
        end
        f_in = [];
    end
    
    T.(nameAroma) = f_aroma;
end

save(fullfile(S.paths.data,sprintf('settings_2_aroma_n%d_%s.mat',size(T,1),datestr(now,'yyyymmdd'))),'S','T');


%% Add motion parameters

doRP = 1;
if doRP
    % -------------------------
    % Add motion parameters 700
    % -------------------------
%     f_rp        = rdir(fullfile(S.paths.rawmri700,'data/aamod_realignunwarp_00001/CBU*/resting/rp_*.txt'));f_rp = {f_rp.name}';
%     f_epiMean   = rdir(fullfile(S.paths.rawmri700,'data/aamod_realignunwarp_00001/CBU*/resting/mean*.nii'));f_epiMean = {f_epiMean.name}';
%     f_rp        = rdir(fullfile(S.paths.rawmri700,'data_fMRI/aamod_realignunwarp_00001/CC*/Rest/rp_*.txt'));f_rp = {f_rp.name}';
%     f_epiMean   = rdir(fullfile(S.paths.rawmri700,'data_fMRI/aamod_realignunwarp_00001/CC*/Rest/mean*.nii'));f_epiMean = {f_epiMean.name}';
    
%     str_cc     = 'CC\d\d\d\d\d\d';
%     SubID       = regexp(f_rp,str_cc,'match','once');
%     Trp         = cell2table([f_rp,f_epiMean,SubID], 'VariableNames',{'f_rp','f_epiMean','SubCCIDc'});
    rel_rms_mean = nan(size(T,1),1);
    parfor isub = 1:size(T,1)
        f_rp = T.f_rpRestW1{isub};
        f_epiMean = regexprep(T.f_restingNormW1{isub},'wr','wart_mean_r');
        
        % Save mean image if it doesnt exist yet
        if ~exist(f_epiMean,'file') 
            V = spm_vol(T.f_restingNormW1{isub});
            Y = spm_read_vols(V);
            Ymean = mean(Y,4);
            Vtemp = V(1);
            Vtemp.fname = f_epiMean;
            spm_write_vol(Vtemp,Ymean);
        end
        [rel_rms1, ~] = y_FD_Jenkinson(f_rp,f_epiMean);
        rel_rms_mean(isub) = mean(rel_rms1);
    end
    T.motion_restingW1 = rel_rms_mean;
    
    % Add WaveletDespike Spike Outlies (outlying wavelet coefficients)
    % to show correlation between spike outliers and movement parameters
    wdsp_mean = nan(size(T,1),1);
    wdsp_median = nan(size(T,1),1);
    parfor isub = 1:size(T,1)
        sp = load(regexprep(T.f_restingWaveletW1{isub},'wds.nii','SP.txt'));
        wdsp_mean(isub) = mean(sp);
        wdsp_median(isub) = median(sp);
    end
    T.motion_restingW1_wdsMean = wdsp_mean;
    T.motion_restingW1_wdsMedian = wdsp_median;    
end


