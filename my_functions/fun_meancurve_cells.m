function [m,e]=fun_meancurve_cells(trials,ch)
    z = nan(length(trials),size(trials{1},2));
    for i=1:length(trials)
        z(i,:) = trials{i}(ch,:);
    end
    m = mean(z);
    e = std(z)./sqrt(size(z,1));
    
end