function [handles] = fun_plotbraintopo_glossy(surf,topo,defaultval,mycfg)

% function plotbraintopo(surf,topo,,defaultval,mycfg)
% surf must be a 3-d brain model including fields '.pnt',and 'tri'
% topo must be a COLUMN vector (ie. extending horizontally, one row many columns)
% defaultval is an arbitrarily to choose from scaling factor

if nargin < 3 
    defaultval = 1;
elseif isempty(defaultval)
    defaultval = 1;
end;
if nargin < 4
    mycfg = struct;
end;
if isfield(mycfg,'mask')
    if isnumeric(mycfg.mask)
        msk = 0;
    elseif strcmp(mycfg.mask,'no')
        msk = 0;
        mycfg = rmfield(mycfg,'mask');
    elseif strcmp(mycfg.mask,'yes')
        msk = 1;
        mycfg = rmfield(mycfg,'mask');
    end;
else
    msk = 0;
end;
        

dumsource.avg.pow = topo;
if (defaultval == 1)
    cfg.colormapping = [min(dumsource.avg.pow),max(dumsource.avg.pow)];
    cfg.opacity = cfg.colormapping;
    if msk == 1
        %cfg.mask = abs(dumsource.avg.pow)';
        cfg.mask = zeros(size(dumsource.avg.pow'));
        cfg.mask(find((dumsource.avg.pow ~=0) & ~isnan(dumsource.avg.pow))) = 1;
    end;
elseif (defaultval ~= 0)
    cfg.colormapping = [-1,1]*defaultval;
    cfg.opacity = cfg.colormapping;
    if msk == 1
        %cfg.mask = abs(dumsource.avg.pow)';
        cfg.mask = zeros(size(dumsource.avg.pow'));
        cfg.mask(find((dumsource.avg.pow ~=0) & ~isnan(dumsource.avg.pow))) = 1;
    end;
elseif (defaultval == 0)
    cfg = struct;
    cfg.colormapping = [-1,1]*max(abs(dumsource.avg.pow));
    cfg.opacity = cfg.colormapping;
end;

if ~isfield(mycfg,'colormapping')
    mycfg.colormapping = cfg.colormapping;
end;
if ~isfield(mycfg,'opacity')
    mycfg.opacity = cfg.opacity;
end;

if msk == 1
    mycfg.mask = cfg.mask;
end;

[handles] = myownsurfplot_a(surf,dumsource,mycfg);
