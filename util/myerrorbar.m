function h = myerrorbar(varargin)
%MYERRORBAR wraps some decent formatting around errorbar function
%   Detailed explanation goes here

% p = MipInputParser;
% p.FunctionName='myerrorbar';
% p.addParameter('color',rgb('blue'),@(x)isnumeric(x)||ischar(x));
% p.parseMagically('caller');
% 
% if ischar(color)
%    color = rgb(color);
% end

% also see formaterrorbar

% sep 2022: at some point i deleted the first argument which was a color
% triplet, but a filesearrch turned up no examples of calling this using
% that syntax, so I added the check below that checks if varargin{1} is a
% color triplet and if so, it should replciate the odl behaviro
if all(size(varargin{1}) == [1 3])
   c  =  varargin{1};
   h  =  errorbar(varargin{2:end});
else
   h  =  errorbar(varargin{:});
   c  =  h.Color;
end

set(h,  ...
      'Marker',           'o',        ...
      'MarkerSize',       8,          ...
      'LineWidth',        1.5,        ...
      'LineStyle',        'none',     ...
      'Color',            c,          ...
      'MarkerFaceColor',  c,          ...
      'MarkerEdgeColor',  'none'      );

end

