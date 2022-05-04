
clear
cd D:\Projects\MSc_AttentionPrediction\MEG\analysisroutines
mysurf = load('canonical_surf_MNI_2052');
mysurf.surf = mysurf.mysurf;
mysurf = rmfield(mysurf,'mysurf');
%gradient:
y = [-100:5:-15];
z = -y*0.7-70;
x = [-90:10:90];
myvoxels = [];
for iy=1:length(y)
    for ix = 1:length(x)
        [~,dumvox,~] = findregionalvoxels(mysurf.surf,[x(ix),y(iy),z(iy)],5,35);
        myvoxels = cat(1,myvoxels,dumvox);
        myvoxels = unique(myvoxels);
    end;
end;
[~,occipitall,~] = findregionalvoxels(mysurf.surf,[-25  -100   10],5,35); %-11.3504  -98.2999   19.8173
[~,occipitalr,~] = findregionalvoxels(mysurf.surf, [25  -100   10],5,35); % 13.1371  -98.2862   16.3996
occipital = union(occipitall,occipitalr);
relevantgridpts = union(myvoxels,occipital);
topo = zeros(size(mysurf.surf.pnt,1),1);
topo(relevantgridpts)=1;
plotbraintopo_glossy(mysurf.surf,topo');
view(100,0);

redundantgridpts = setdiff([1:size(mysurf.surf.pnt,1)],relevantgridpts);
reducedgrid = mysurf.surf;
reducedgrid.pnt(redundantgridpts,:) = 0;

save('ventral_occ_tempgrid_2052','reducedgrid','relevantgridpts','redundantgridpts');



    