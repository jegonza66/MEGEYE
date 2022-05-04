function [m,e]=fun_meancurve(timelock,indch)
    if ~isfield(timelock,'trial')
        display('ERROR: Field trial is left. Calculate the timelock again and check the cfg.')
        return
    end
    z = squeeze(mean(timelock.trial(:,indch,:),2));
    m = mean(z);
    e = std(z)./sqrt(size(z,1));
    
end
    
    