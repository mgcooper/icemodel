function [handle] = addOnetoOne(varargin )
%addOnetoOne adds a one:one line to an open plot
%   Adds a solid black line of width 1, scales the x and y lims to be
%   equal, returns the handle to the plot so the user can specify line
%   properties on the 1:1 line. Alternatively, the line properties can be
%   passed in i.e. "varargin"

% INPUTS: 
% optional, can be left blank, or any standard line property can be used as
% an input to modify the 1:1 line. See examples.

% OUTPUTS: 
% handle = handle to the 1:1 line object

% AUTHOR: Matthew Guy Cooper, UCLA, guycooper@ucla.edu

% EXAMPLE 1: simplest case, just add the one to one line
% x = 1:1:100;
% y = x.*(rand(1,100) + .5);
% h = scatter(x,y);
% addOnetoOne;

% EXAMPLE 2: change the default 1:1 line to a black dashed line
% x = 1:1:100;
% y = x.*(rand(1,100) + .5);
% h = scatter(x,y);
% addOnetoOne('k--')

% EXAMPLE 3: change the default 1:1 line to very thick line;
% x = 1:1:100;
% y = x.*(rand(1,100) + .5);
% h = scatter(x,y);
% addOnetoOne('k--','linewidth',5)

% if nargin==0
%     posneg = pos;
% elseif nargin>0
%     if strcmp(posneg,'positive')
%     end
% end

hold on; axis tight;

xlims       =   xlim;
ylims       =   ylim;

newlims(1)  =   min(min(xlims),min(ylims));
newlims(2)  =   max(max(xlims),max(ylims));

set(gca,'XLim',newlims,'YLim',newlims);

delta       =   (newlims(2) - newlims(1))/100;

handle      =   plot(newlims(1):delta:newlims(2),newlims(1):delta:newlims(2), ...
                        varargin{:});

% i diabled this b/c it's the reason the exponent isn't showing up on the
% tick labels anymore, need to come up with a solution 
% % now update the ticks and labels
% xticks      =   get(gca,'Xtick');
% yticks      =   get(gca,'Ytick');
% if length(xticks) >= length(yticks)
%     ticks       =   xticks;
%     ticklabels  =   get(gca,'XTickLabel');
% else
%     ticks       =   yticks;
%     ticklabels  =   get(gca,'YTickLabel');
% end
% 
% set(gca,'XTick',ticks,'XTickLabel',ticklabels,'YTick',ticks, ...
%         'YTickLabel',ticks);

axis square
end


% this is the original version:
% hold on;
% 
% xlims       =       get(gca,'XLim');
% ylims       =       get(gca,'YLim');
% newlims(1)  =       min(xlims(1),ylims(1));
% newlims(2)  =       max(xlims(2),ylims(2));
%                     set(gca,'XLim',newlims,'YLim',newlims);
% handle      =       plot(newlims(1):.1:newlims(2),newlims(1):.1:newlims(2),varargin{:});
% axis square





































