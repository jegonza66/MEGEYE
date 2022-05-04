
function [myfilters,filterselect] = findbeamfcoeffs4ROIlowreshiresgrid(mysurflo,mysurfhi,source,searchregion)

% [myfilters,filterselect] = findbeamfcoeffs4ROIlowreshiresgrid(mysurflo,mysurfhi,source,searchregion)
%
% - myfilters:  gives a cell array for the filters of all gridpoints specified
%               in mysurflo, as they can be found in mysurfhi within the searchregion
% - filterselect:   gives the indices of those filters for the exactly
%                   corresponding gridpoints in mysurfhi to those of mysurflo

    myfilters = [];
    filterselect = zeros(size(mysurflo.pnt,1),1);
    for igrid = 1:size(mysurflo.pnt,1)
        [neighbrgidpts,dum] = findregionalvoxels(mysurfhi,mysurflo.pnt(igrid,:),searchregion,0.1);
        filterselect(igrid,1) = dum(1);
        for k=1:length(neighbrgidpts)
            myfilters{igrid}(k,:) = source.avg.filter{neighbrgidpts(k)}(1,:);
        end;

    end;

end
