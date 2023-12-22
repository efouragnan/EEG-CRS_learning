function [ ] = setFigDefaults( )
% gives nice default figure properties
%
% check default by:        get(groot, 'default')
% check default fig pos:   get(0,'defaultfigureposition')
% get back to defaults:    set(groot, 'Default', struct())
%

defMark        = 5;
defaultsize    = 11;     % everything else
smallfontsize  = 9;
defLinewidth   = 1.5;
fontmulti      = 1.1;



% general properties of figure;:
set(groot,'defaultAxesbox','off'); % does not work because overwritten by plot command
set(groot,'defaultFigureColor','w');
set(groot,'defaultFigureclipping','off');

% default text:
set(groot,'defaultAxesfontsize',defaultsize);
%set(groot,'defaultAxesLinewidth',defLinewidth);
set(groot,'defaultLineLinewidth',defLinewidth);
set(groot,'defaultLineMarkerSize',defMark)
set(groot,'defaultAxesclipping','off');
set(groot,'defaultAxesLayer','top');
set(groot,'defaultAxesFontName','Arial') ; %,'PlotBoxAspectRatio',[4/3,1,1])
set(groot,'defaultAxesTickDir','out');
%set(groot,'defaultAxesTickLength',[1 1]*0.02); % leave it to default right now


% title and axes
set(groot,'defaultAxeslabelFontsizemultiplier',fontmulti);  % for some reason Fontsize is intepreted as Fontsizemultiplier
set(groot,'defaultAxesTitleFontsizemultiplier',fontmulti);
set(groot,'defaultTextinterpreter','none');
set(groot,'defaultAxesXColor','k');                                               % the color of the axis line and the tick marks
set(groot,'defaultAxesYColor','k');
set(groot,'defaultAxesticklabelinterpreter','none'); 

% adjust legend:
set(groot,'DefaultLegendfontsize',smallfontsize); 
set(groot,'DefaultLegendbox','on');
set(groot,'DefaultLegendinterpreter','none');
set(groot,'DefaultAxesTickDir','out');

%mpos        = get(groot,'monitorpositions'); if size(mpos,1)>1, mpos=mpos(1,:); end;
%set(groot,'defaultfigureposition',mpos); % make screen size the default

end

