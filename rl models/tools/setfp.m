function [ ] = setfp(fighandle)
%setfp(gcf)
% sets figure properties, only needed once per figure;  short version of set_default_figure_properties, bc everything else is defined in setFigDefaults configuring figure property defaults;

% this is only necessary bc "plot" overwrites defaults.

allaxes = findall(fighandle,'type','axes');

for ia=1:numel(allaxes)
   axishandle=allaxes(ia);
   set(axishandle,'tickdir','out'); 
   set(axishandle,'box','off');
   set(axishandle,'Clipping','on');
end


end