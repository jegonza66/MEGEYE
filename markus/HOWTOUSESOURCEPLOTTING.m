mysurf = load('canonical_surf_MNI_413');
mysurf.surf.pos = mysurf.surf.pnt;
gridd.pos = mysurf.surf.pnt;
centralsensmot = find((gridd.pos(:,2) < 15) & (gridd.pos(:,2) > -60) & (gridd.pos(:,3) > 25));
topo = zeros(413,1);
topo(centralsensmot)=1;
pcfg.material = 'shiny';
pcfg.lightangle = [6,-16]; % [-8,-6];
plotbraintopo_glossy(mysurf.surf,topo',[10],pcfg);