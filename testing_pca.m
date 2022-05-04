load(cfg.datafname,'data');
cocfg = [];
cocfg.blc           = 'yes';
cocfg.method        = 'pca';
cocfg.channel       = {'all'};
cocfg.numcomponent  = length(data.label);
comp                = ft_componentanalysis(cocfg,data);

 
%   Instead of specifying a component analysis method, you can also specify
%   a previously computed unmixing matrix, which will be used to estimate the
%   component timecourses in this data. This requires
%     cfg.unmixing     = NxN unmixing matrix
%     cfg.topolabel    = Nx1 cell-array with the channel labels

