function [S,T] = funanat_03_parcellate(S,T,dataname,wave,parcellation)
% Function to parcellate pre-processed fMRI time-series (AROMA outputs)


nameParc = parcellation.name;

Ns = size(T,1);

% ---------------------------
% Predefining some variables
% ---------------------------
varname = ['resting' dataname 'W' num2str(wave)];
f_in_signal = T.(['f_' varname]);
f_in_roi    = T.(['f_' varname]);


outdir  = sprintf('%s/restingW%d/data_node_%s/',S.paths.postproc,wave,lower(dataname));
mkdir(outdir);

switch nameParc
      case 'hcp360'%  HCP Glasser Parcellation
        f_ROImask = parcellation.f_mask;
        V = spm_vol(f_ROImask);
        
        [ROImask BNinfo] = ParseInNii(f_ROImask,'compress',0);% load GM Dartel mask
        ROIlabel  = S.parcellation.hcp360.info;
        ROIlabel.ROIname = ROIlabel.region_short_name;
        ROIlabel.ROIval  = ROIlabel.region_no;
        Nroi             = size(ROIlabel,1);
end

% ---------------------------------------------------------
% Some settings for extracting globals signals
% ---------------------------------------------------------
csfMask     = spm_read_vols(spm_vol(S.process.f_csfMask))   > 0.8;
wmMask      = spm_read_vols(spm_vol(S.process.f_wmMask))    > 0.8;
gmMask      = spm_read_vols(spm_vol(S.process.f_gmMask))    > 0.8;
brainMask   = spm_read_vols(spm_vol(S.process.f_brainMask)) > 0.8;

Signals = {'csf'    csfMask     0.8;...
           'wm'     wmMask      0.8;...
           'gm'     gmMask      0.8;...
           'brain'  brainMask   0.8;...
           };

%%
% Loop through all subjects
% --------------------------
parfor isub = 1:Ns

    if ~isempty(f_in_roi{isub})

    f_out = fullfile(outdir,T.SubID{isub},sprintf('node_%s_%s_%s.mat',nameParc,T.SubID{isub},varname));
    
        if ~exist(f_out,'file')
            mkdir(fileparts(f_out));
            B = ParseInNii(f_in_roi{isub},'compress',0);
            ROI = struct();
            
            % Loop through all ROIs
            for iroi = 1:Nroi
                roi_val     = ROIlabel.ROIval(iroi);
                idx_vox     = ROImask==roi_val;
                dat         = B(:,idx_vox);
                meanSignal  = mean(dat,2);
                [U,latent]  = svd(dat);
                explained   = 100*latent/sum(latent);%variances of all individual principal components

                % -------------------------------------------------------------
                % Ensure 1st PC timecourse is positively correlated with mean
                % timecourse (SVD/PCA can flip the sign)
                % -------------------------------------------------------------
                r           = corr(U(:,1),meanSignal);
                pcaSignal   = U(:,1)*r;
                ROI(iroi).name       = ROIlabel.ROIname{iroi};
                ROI(iroi).id         = ROIlabel.ROIval(iroi);
                ROI(iroi).meanSignal = meanSignal;
                ROI(iroi).pcaSignal  = pcaSignal;
                ROI(iroi).explained  = explained(1:5);
            end
            
           
            B = ParseInNii(f_in_signal{isub},'compress',0);
            
            Tissue = struct();
            for isignal = 1:size(Signals,1)
                nameSignal = Signals{isignal,1};
                dat        = B(:,Signals{isignal,2}(:));
                meanSignal =  mean(dat,2);
                [U,latent] = svd(dat);
                explained  = 100*latent/sum(latent);%variances of all individual principal components
                pcaSignal  =  U(:,1:5);
                Tissue.(nameSignal).meanSignal = meanSignal;
                Tissue.(nameSignal).pcaSignal  = pcaSignal;
                Tissue.(nameSignal).explained  = explained(1:5);
            end
            
            kat_parfor_save(f_out,{'ROI','Tissue'},ROI,Tissue);

            
        end
        f_out_all{isub} = f_out;
    end
end

T.(['f_node_' varname '_' nameParc]) = f_out_all';