function [withinsphere,roiindx,withincube] = findregionalvoxels(grid,roi,dist,roidist)

% function [withinsphere,roiindx,withincube] = findregionalvoxels(grid,roi,dist,roidist)
%
% grid is a structure containing a position structure
% roi is a row vector defining the center of the ROI
% dist is the search radius
%
% withinsphere,withincube are indices to voxels inside the cube / sphere of defined 'radius'

if nargin < 4
    roidist = 2;
end;

if ~isfield(grid,'pos')
    grid.pos = grid.pnt;
end;

surr_x = find(abs(grid.pos(:,1)-roi(:,1)) < dist);
surr_y = find(abs(grid.pos(:,2)-roi(:,2)) < dist);
surr_z = find(abs(grid.pos(:,3)-roi(:,3)) < dist);

withincube = intersect(surr_x,surr_y); 
withincube = intersect(withincube,surr_z);

distance = sqrt(sum((grid.pos(withincube,:) - repmat(roi(1,:),[length(withincube),1])).^2,2));
withinsphere = find(distance < dist);
withinsphere = withincube(withinsphere);

surr_x = find(abs(grid.pos(:,1)-roi(:,1)) < roidist);
surr_y = find(abs(grid.pos(:,2)-roi(:,2)) < roidist);
surr_z = find(abs(grid.pos(:,3)-roi(:,3)) < roidist);

roiindx = intersect(surr_x,surr_y); 
roiindx = intersect(roiindx,surr_z);
