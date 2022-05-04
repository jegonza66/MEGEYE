
function [handle] = myownsurfplot(surf,source,cfg)

% function [handles] = myownsurfplot(surf,source,cfg)
%
% allows to plot source results obtained from a cortical mesh
%
% cfg - opacity: specifies the opacities of the vertices e.g. [0,1]
% cfg - colormapping: defines the colorrange of the vertices e.g. [0,1]
% cfg - edgealpha: defines the transparency of the edges for identifying gyri & sulci, e.g. 0.2
% cfg - mask matrix specifying opacity, format: [n,1] for n vertices

if nargin < 3
    cfg = struct();
end;

    
% cortex_light = [0.781 0.762 0.664];
% cortex_dark  = [0.781 0.762 0.664]/2;
cortex_light = [1 1 1]*0.9;
cortex_dark  = [1 1 1]*0.9/1.5;

  
sourcevals = source.avg.pow(:);
backgcolor = repmat(cortex_light, size(surf.pnt,1), 1);

if ~isfield(cfg,'opacity')
    opacmin = min((source.avg.pow(:)));
    opacmax = max((source.avg.pow(:)));
else
    opacmin = cfg.opacity(1);
    opacmax = cfg.opacity(2);
end;
if ~isfield(cfg,'colormapping')
    fcolmin = min((source.avg.pow(:)));
    fcolmax = max((source.avg.pow(:)));
else
    fcolmin = cfg.colormapping(1);
    fcolmax = cfg.colormapping(2);
end;
if ~isfield(cfg,'edgealpha'), edgealpha = 1; end;
if ~isfield(cfg,'mask'), 
    maskval = ones(size(source.avg.pow(:)));
    hasmsk = 0;
else
    maskval = cfg.mask;
    hasmsk = 1;
end;

if ~isfield(cfg,'handle')
    handle = figure;
end
h1 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', backgcolor , 'FaceColor', 'interp');
%this has changed
set(h1, 'EdgeColor', 'none');
%set(h1, 'EdgeColor', [0,0,0],'EdgeAlpha',edgealpha);
axis   off;
axis vis3d;
axis equal;

h2 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', sourcevals , 'FaceColor', 'interp');
set(h2, 'EdgeColor',  'none');
%set(h2, 'EdgeColor',(cortex_light+cortex_dark)/2,'EdgeAlpha',0.15);%

alphamap('vdown');

if hasmsk
%     set(h2, 'AlphaData', maskval);
    alim('manual');
    %set(h2, 'AlphaDataMapping','scaled');
    set(h2, 'AlphaDataMapping','none');
    set(h2, 'FaceAlpha',   'interp');
    set(h2, 'FaceVertexAlphaData', maskval);
    %alim(gca, [opacmin opacmax]);
end
caxis(gca,[fcolmin fcolmax]);
%alim(gca, [opacmin opacmax]);


if isfield(cfg,'material')
    if strcmp(cfg.material,'shiny')
        lighting phong;
        material shiny;
        
    else
        lighting gouraud; %phong %
        material dull ;
    end;
else
        lighting gouraud; %phong %gouraud;
        material dull ;
end;    
        

if ~isfield(cfg,'lightangle')
    camlight %('style','infinite')  %light 
else
    lightangle(cfg.lightangle(1),cfg.lightangle(2));
end;

colormap('jet');
colorbar;

set(gcf,'Color','w')

% handles.h1 = h1;
% handles.h2 = h2;


