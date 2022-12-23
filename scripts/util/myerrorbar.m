function h = myerrorbar(varargin)
%MYERRORBAR wraps some decent formatting around errorbar function
% 
%   MYERRORBAR(Y,E) plots Y and draws a vertical error bar at each element of
%   Y. The error bar is a distance of E(i) above and below the curve so
%   that each bar is symmetric and 2*E(i) long.
%
%   MYERRORBAR(X,Y,E) plots Y versus X with symmetric vertical error bars
%   2*E(i) long. X, Y, E must be the same size. When they are vectors, each
%   error bar is a distance of E(i) above and below the point defined by
%   (X(i),Y(i)). When they are matrices, each error bar is a distance of
%   E(i,j) above and below the point defined by (X(i,j),Y(i,j)).
%
%   MYERRORBAR(X,Y,NEG,POS) plots X versus Y with vertical error bars
%   NEG(i)+POS(i) long specifying the lower and upper error bars. X and Y
%   must be the same size. NEG and POS must be the same size as Y or empty.
%   When they are vectors, each error bar is a distance of NEG(i) below and
%   POS(i) above the point defined by (X(i),Y(i)). When they are matrices,
%   each error bar is a distance of NEG(i,j) below and POS(i,j) above the
%   point defined by (X(i,j),Y(i,j)). When they are empty the error bar is
%   not drawn.
%
%   MYERRORBAR( ___ ,Orientation) specifies the orientation of the error
%   bars. Orientation can be 'horizontal', 'vertical', or 'both'. When the
%   orientation is omitted the default is 'vertical'.
%
%   MYERRORBAR(X,Y,YNEG,YPOS,XNEG,XPOS) plots X versus Y with vertical error
%   bars YNEG(i)+YPOS(i) long specifying the lower and upper error bars and
%   horizontal error bars XNEG(i)+XPOS(i) long specifying the left and
%   right error bars. X and Y must be the same size. YNEG, YPOS, XNEG, and
%   XPOS must be the same size as Y or empty. When they are empty the error
%   bar is not drawn.
%
%   MYERRORBAR( ___ ,LineSpec) specifies the color, line style, and marker.
%   The color is applied to the data line and error bars. The line style
%   and marker are applied to the data line only.
%
%   MYERRORBAR(AX, ___ ) plots into the axes specified by AX instead of the
%   current axes.
% 
%   MYERRORBAR(C, ___ ) uses the rgb color triplet C instead of the default
%   colororder. (mgc)
%
%   H = MYERRORBAR( ___ ) returns handles to the errorbarseries objects
%   created. ERRORBAR creates one object for vector input arguments and one
%   object per column for matrix input arguments.
%
%   Example: Draws symmetric error bars of unit standard deviation.
%      x = 1:10;
%      y = sin(x);
%      e = std(y)*ones(size(x));
%      errorbar(x,y,e)

% p = magicParser;
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