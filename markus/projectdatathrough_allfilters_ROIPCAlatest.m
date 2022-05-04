function [newdat] = projectdatathrough_allfilters_ROIPCAlatest(data,config)

filter = config.filter;
%assuming the filters for one ROI are in a cell array, with individual grid
%points in rows and channel weights in columns
% this function assumes the filters coming from
% "findbeamfcoeffs4ROIlowreshiresgrid.m"

%do make sure the channels of source structure (filter) and data
%correspond!!!
if isfield(config,'sourcelabel')  
    [datasel,filtersel] = match_str(data.label,config.sourcelabel);
else
    datasel = [1:length(data.label)];
    filtersel = [1:size(myfilters{1},2)];
end;

%create an artificial dataset for obtaining the PC coefficients for each
%filter/ROI for the entire dataset - use only a subset of trials to make
%PCA faster
trlselection = [];
if length(data.trial) > 400
    trlsel = [1:2:length(data.trial)];
elseif length(data.trial) > 200
    trlsel = [1:1:length(data.trial)];  
else
    trlsel = [1:length(data.trial)];      
end;
selcfg.trials = trlsel;
mydata = ft_selectdata_new(selcfg,data);

%before PCA - filter the helper data according to signal of interest
if isfield(config,'preproc')
    ppcfg = config.preproc;
    mydata = ft_preprocessing(ppcfg,mydata);
end;

%obtain the PC coefficients for the entire dataset - for each ROI
datdum = cell2mat(mydata.trial)'; %has channels now in columns
clear mydata;
filterpcacoeffs = {};
for z=1:length(filter)
    %for each low res grid point will now have many filters
    % reduce to prevent pca 
    if size(filter{z},1) > 40
        filter{z} = filter{z}([1:4:end],:);
    end;
    dumdat = datdum(:,datasel) * filter{z}(:,filtersel)';
    %dumdat = nanzscore(dumdat,1);
    coeff = princomp(dumdat);
    filterpcacoeffs{z} = coeff(:,1);
end;
clear datdum dumdat;

%finally create the virtual channel timeseries
newdat = rmfield(data,{'trial';'label';'grad'});
dumdat = [];
for trl = 1:length(data.trial)
    for z=1:length(filter)
        %for each low res grid point will now have many filters
        dumdat = (data.trial{trl}(datasel,:)' * filter{z}(:,filtersel)');
        newdat.trial{trl}(z,:) = dumdat*squeeze(filterpcacoeffs{z});
    end
end;
clear dumdat filterpcacoeffs;

for q = 1:length(filter)
    label{q,1} = ['virtfilt_',num2str(q)];
end;
newdat.label = label;
clear data;

end
