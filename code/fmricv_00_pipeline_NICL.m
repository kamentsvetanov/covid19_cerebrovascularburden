clear

% ----------------------
% Dependencies 
% -----------------------
% fdr_bky - https://github.com/kamentsvetanov/external/tree/master/mat/matlabcentral

% Set some paths
S.paths.figures = '/home/kt03/Projects/public-code/covid19_cerebrovascularburden/results';

% ----------------------------------------------------------
% 2. Pre-process fMRI timeseries using AROMA/Wavelet
% ----------------------------------------------------------
[S,T] = fmricv_02_process_resting_NICL(S,T);

% ----------------------------------------------------------
% 3. Parcellate pre-processed fMRI timeseries
% ----------------------------------------------------------
[S,T] = fmricv_03_parcellate_NICL(S,T,'Aroma',1,S.parcellation.hcp360);

% -------------------------------------------------------------------------
% 4. Post-process motion-corrected 4D images and compute RSFA
% -------------------------------------------------------------------------
[S,T] = fmricv_04_postproc_func_4D_NICL(S,T,'Aroma',1,'pcaSignal',1);

%% Identify CV19-related IC of RSFA
T=readtable('/home/kt03/Projects/public-code/covid19_cerebrovascularburden/data/data_patients_controls.xlsx');

% Determine effects of interest across all components
% ---------------------------------------------------
DV    = 'loadings_rsfa'; % Name of depenent variable
IV    = 'Age';      % Independent variable
MV    = 'GroupIdx'; % Moderator
CV    = {'Sex'};    % Covariates of no interest
numIC = 10;         % Number of ICs
SS    = vartype('numeric');   

for icomp = 1:numIC
    tbl = table();
    tbl = T(:,[IV MV CV]);
    DVi = [DV '_' num2str(icomp)]; 
    DVname = ['IC' num2str(icomp)];
    tbl.(DVname) = T.(DVi);
    tbl(:,SS) = normalize(tbl(:,SS));
    
    model = sprintf('%s ~ %s*%s + %s', DVname,IV,MV, strjoin(CV, ' + '));
    mlr = fitlm(tbl,model,'RobustOpts',1);
    MLR{icomp} = mlr;
    pvals(:,icomp) = mlr.Coefficients.pValue; 
end

id_coeff = 3; % Which coefficient is Group effect
id_comp  = find(fdr_bky(pvals(3,:))); % Components showing significant groupp differences, FDR-corrected
MLR{id_comp}
fprintf('\n \n IC %d shows significant group differences (fdr-corrrected). \n',id_comp);

%% Bar plot of IC subject loadings
icomp       = 4;

% Remove the effects of age to show to visualise correctly the group
% effects on RSFA as in the linear model above.
dat = T.([DV '_' num2str(icomp)]);
mlr = fitlm([T.Age ],dat);

cfg.dat                 = zscore(mlr.Residuals.Standardized);
cfg.group_idx           = T.GroupIdx+1;
cfg.doFlipRelationship  = 1; % Flip Loadings so that higher loading means higher CVR burden
cfg.doSaveFigure        = 0;
cfg.group_names         = {'Controls';'Patients'};
cfg.condition_names     = sprintf('IC%03d',icomp);
cfg.fname               = fullfile(S.paths.figures,sprintf('Figure2_boxplot_groupEffect_%s_IC%03d.png',DV,icomp));
fmricv_util_plotting_boxplot_groupEffects_NICL(cfg);

%% 	Linking COVID-19 severity to cerebrovascular impairment and their correlates
clear;
Tpatient = readtable('/home/kt03/Projects/public-code/covid19_cerebrovascularburden/data/data_patients.xlsx');
SS   = vartype('numeric');


% Correlation between COVID-19 Severity and cerebrovascular impairment
% using PLS saliences from first significant component
% -------------------------------------------------------------------------
fprintf('Correlation between saliences is %.3f. \n ',corr(Tpatient.RSFA,Tpatient.CV19Severity));


% Test if other factors explain the relationship between RSFA and CV19
% -------------------------------------------------------------------------
model = 'RSFA ~ Sex + Age + Clinical_PCA1_iAll + Clinical_PCA2_iAll + CV19Severity' ;
mlr = fitlm(normalize(Tpatient(:,SS)), model,'RobustOpts',1)

% Test if comorbidity factors influence RSFA-CV19 relationship 
% -------------------------------------------------------------------------
modelComorbidity = [model '+ comorbidity_cardiovascular+comorbidity_diabetes+comorbidity_hypertension+comorbidity_respiratory'];
mlrComorbiity = fitlm(normalize(Tpatient(:,SS)), modelComorbidity,'RobustOpts',1)


% Identify commonality between other factors and RSFA-CV19 link
% Requires commonality analysis for neuroimaing toobox
% https://github.com/kamentsvetanov/CommonalityAnalysis
% -------------------------------------------------------------------------
cfg         = [];
cfg.mlr     = mlr;
cfg.doPerm  = 1;
cfg.numPerm = 10000;
cfg.runParfor = 1;
CA1 = ca_stats_commonality(cfg)

% COVID-19-related cerebrovascular burden association with physical, 
% cognitive, and mental functioning (PCM PC1 and PCM PC2)
% -------------------------------------------------------------------------

% No evidence that RSFA is related to PCM 1
model = ['Cogment_PCA1_iAll ~ Age+ Sex + RSFA']; 
mlr = fitlm(normalize(Tpatient(:,SS)), model,'RobustOpts',1)

% RSFA is significantly related to PCM1
model = ['Cogment_PCA2_iAll ~ Age+ Sex + RSFA'];
mlr = fitlm(normalize(Tpatient(:,SS)), model,'RobustOpts',1)

% RSFA is significantly related to PCM1 over and above commorbidity factors
modelComorbidity = [model '+ comorbidity_cardiovascular+comorbidity_diabetes+comorbidity_hypertension+comorbidity_respiratory'];
mlrComorbidity = fitlm(normalize(Tpatient(:,SS)), modelComorbidity,'RobustOpts',1)

%% Spatial overlap of cv19severity with receptor, metabolic statmaps and Covid19-related genes - Figure 4
% -----------------------------------------------------    
clear
load('/home/kt03/Projects/public-code/covid19_cerebrovascularburden/data/results_spatial_overlap.mat');

% -------------------------------------------------------------------------
% Add target class to the table consitent with Hansen et al 2021, SI F2
mapclass    = {'serotonin',     {'5HT'};...
               'domanine',      {'D1','D2','DAT'};...
               'norepinephrine',{'NAT'};...
               'histamine',     {'H3'};...
               'acetylcholine', {'A4B2','M1','VAChT'};...
               'cannabinoid',   {'CB1'};...
               'opioid',        {'MU'};...
               'glutamate',     {'mGluR5'};...
               'GABA',          {'GABAa-bz'};...
               'metabolism',    {'GI','CMRO2','CMRGlu','CBV','CBF',};...
               'gene expression',{'CTSB','CTSL','IFNAR2','LY6E','NRP2','NRP1','ACE2'};...
               };
Mcv19.class = strings(height(Mcv19),1);
Mcv19.class_sort = nan(height(Mcv19),1);
for iclass = 1:size(mapclass)
    idx = contains(Mcv19.target,mapclass{iclass,2});
    Mcv19.class(idx) = mapclass{iclass,1};
    Mcv19.class_sort(idx) = iclass; 
end
[Mcv19,I] = sortrows(Mcv19,{'class_sort','target'});


% -------------------------------------------------------------------------
% Minor relabelling of target and assign as labels
labels = Mcv19.target;
labels = regexprep(labels,'5HT','5-HT');
labels = regexprep(labels,'NAT','NET');
labels = regexprep(labels,'MU','MOR');
labels = regexprep(labels,'GABAa-bz','GABAa/bz');
labels = regexprep(labels,'CV19-gene-','');
Mcv19.label = labels;

rvals     = Mcv19.rmap;
ncolor    = 256;
cmap      = cbrewer('seq','YlOrRd',ncolor);
facecolor = cmap(ceil(ncolor*.4),:);
facecolor = [255 214 165]./255;
edgecolor = [224,122,0]./255;
figure('Position',[100 100 300 1600]);
for ibar = height(Mcv19):-1:1
    bh(ibar) = barh(ibar,abs(rvals(ibar)),0.4);
    if rvals(ibar)<0
        linestyle = '--';
    else
        linestyle = 'none';
    end
    set(bh(ibar),'FaceColor',facecolor,'EdgeColor',(edgecolor),'LineStyle',linestyle);% 
    text(0.002,ibar-.4,Mcv19.label{ibar},'FontWeight','light');
    hold on;
end
set(gca,'Ydir','reverse');
set(gca,'Ytick',[]);
set(gca,'box','off','ycolor','w');
set(gca,'Xtick',[0:.2:.6])
set(gca, 'color', 'none');
set(gca,'fontweight','light');
xlabel('spatial correlation');
title('CV19-Severity induced CVB');


