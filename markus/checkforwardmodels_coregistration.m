function [handles] = checkforwardmodels_coregistration(spm_inversestructure)

% function checkforwardmodels_coregistration(spm_inversestructure)
%
% use as checkforwardmodels_coregistration(D.inv{1})
%
% markus bauer

%% PLOT THE MRI from Forward model
[mysurf] = convert_meshintosurf(spm_inversestructure.forward.mesh);
plotcfg.opacity = [-1,1];
plotcfg.colormapping = [0,2];
plotbraintopo_glossy(mysurf,ones(8196,1),2,plotcfg);
set(gcf,'Name','FORWARD mesh');

%plot the volume conductor shell on top of it
hold on;
mypoints = spm_inversestructure.forward.vol.bnd.pnt;
plot3(mypoints(1:3:end,1),mypoints(1:3:end,2),mypoints(1:3:end,3),'.r');
%plot the sensors on top of it
mypoints = spm_inversestructure.datareg.sensors.chanpos;
plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.k');

%plot the fiducials
mypoints = spm_inversestructure.datareg.fid_mri.fid.pnt;
plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'ob','LineWidth',3);
% these are the fiducial positions as originally measured
mypoints = spm_inversestructure.datareg.fid_eeg.fid.pnt;
plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'om','LineWidth',3);
legend('','','','','MRIbased fid','MEG measured fid');

handles(1) = gcf;


% %% plot the surface of boundary element with fiducials
% [mysurf] = convert_meshintosurf(spm_inversestructure.forward.mesh);
% plotcfg.opacity = [-1,1];
% plotcfg.colormapping = [0,2];
% plotbraintopo_glossy(mysurf,ones(8196,1),2,plotcfg);
% set(gcf,'Name','FORWARD shell');
% %plot the sensors on top of it
% hold on;
% mypoints = spm_inversestructure.datareg.sensors.chanpos;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.k');
% %plot the fiducials
% mypoints = spm_inversestructure.datareg.fid_mri.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'ob');
% % these are the fiducial positions as originally measured
% mypoints = spm_inversestructure.datareg.fid_eeg.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'or');
% legend('MRIbased fid','MEG measured fid');
% 
% handles(3) = gcf;
% 
% %% PLOT THE TESSELATED MRI - standard brain or just different coordinates??
% [mysurf] = convert_meshintosurf(spm_inversestructure.mesh.tess_mni);
% plotcfg.opacity = [-1,1];
% plotcfg.colormapping = [0,2];
% plotbraintopo_glossy(mysurf,ones(8196,1),2,plotcfg);
% hold on;
% set(gcf,'Name','MNI mesh');
% mypoints = spm_inversestructure.mesh.fid.pnt;
% plot3(mypoints(1:3:end,1),mypoints(1:3:end,2),mypoints(1:3:end,3),'.r');
% mypoints = spm_inversestructure.mesh.fid.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'ob');
% 
% handles(1) = gcf;




% 
% % DO THE SAME FOR "FID_MRI" - assuming this is from the original MNI brain
% mypoints = spm_inversestructure.datareg.fid_mri.fid.pnt;
% figure;plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'ob');
% hold on;
% % these are the fiducial positions as originally measured
% mypoints = spm_inversestructure.datareg.fid_eeg.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'or');
% % these are the MEG sensors
% mypoints = spm_inversestructure.datareg.sensors.chanpos;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.k');
% legend('MR fiducials','MEG fiducials','sensors');
% set(gcf,'Name','fiducials and sensors');
% 
% % PLOT THE head SURFACE used in FORWARD MODEL
% % these are the fiducial positions as used for source reconstruction
% figure
% mypoints = spm_inversestructure.forward.vol.bnd.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.b');
% hold on;
% % these are the fiducial positions as originally measured
% mypoints = spm_inversestructure.datareg.fid_eeg.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'or');
% mypoints = spm_inversestructure.datareg.fid_mri.fid.pnt;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'og');
% % these are the MEG sensors
% mypoints = spm_inversestructure.datareg.sensors.chanpos;
% plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.k');
% set(gcf,'Name','FORWARD MODEL headsurface');
% 
% 
% 
% % 
% % % PLOT THE SURFACE of warped HEADSHAPE WITH its FIDUCIALS
% % mypoints = spm_inversestructure.mesh.fid.pnt;
% % figure;plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.b');
% % hold on;
% % % these are the fiducial positions as used for source reconstruction
% % mypoints = spm_inversestructure.mesh.fid.fid.pnt;
% % plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'og');
% % % these are the fiducial positions as originally measured
% % mypoints = spm_inversestructure.datareg.fid_eeg.fid.pnt;
% % plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'or');
% % % these are the MEG sensors
% % mypoints = spm_inversestructure.datareg.sensors.chanpos;
% % plot3(mypoints(:,1),mypoints(:,2),mypoints(:,3),'.y');
% % set(gcf,'Name','individual headsurface and fiducials');
% % 
% 
