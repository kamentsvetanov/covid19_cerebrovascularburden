function fmricv_util_plotting_boxplot_groupEffects(cfg)
%% Boxplot group differences of ICA loadings

% Unpack cfg
dat                 = cfg.dat;
doFlipRelationship  = cfg.doFlipRelationship;
group_idx           = cfg.group_idx;
group_names         = cfg.group_names;
condition_names     = cfg.condition_names;
doSaveFigure        = cfg.doSaveFigure;
fname               = cfg.fname;


% an alternative color scheme for some plots
c =  [0.45, 0.80, 0.69;...
      0.98, 0.40, 0.35;...
      ];  
 
% figure;
% H = plotAdded(mlr,find(strcmp(mlr.CoefficientNames,'GroupIdx')));
% dat =  H(1).YData;
if doFlipRelationship
   dat = -dat;
   fname = regexprep(fname,'.png','flipped.png');
end

% T.GroupCat = categorical(T.GroupIdx,[1, 0],{'Patients', 'Controls'});

figure;
h = daboxplot(dat,'groups',group_idx,'xtlabels', condition_names,...
    'colors',c,'fill',0,'whiskers',0,'scatter',2,...
    'outliers',0,'scattersize',30,'flipcolors',1,'boxspacing',1,...
    'legend',group_names); 
ylabel('Subject Loadings');
ylim([-3 3]);
% xl = xlim; xlim([xl(1), xl(2)+0.75]); % make more space for the legend
set(gca,'FontSize',12,'Color','none','YTick',[-3:3]);


if doSaveFigure
    export_fig(fname,'-transparent','-painters','-q200','-r600');
end
