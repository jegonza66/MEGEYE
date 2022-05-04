function data = fun_include_eyech(datain,datain_weyech)
%     datain = cfg.datafname;
%     datain_weyech = cfg.datafname_weyech;

    load(datain_weyech,'data')
    datatmp = data;
    clear data
    indch = find(cellfun(@(x) ~isempty(strfind(x,'UADC')),datatmp.label));
    % datatmp.label(indch);
    load(datain,'data')

    for tr=1:length(datatmp.trial)
        data.trial{tr} = [data.trial{tr}; datatmp.trial{tr}(indch,:)];
    end
    data.label = [data.label; datatmp.label(indch)];
end