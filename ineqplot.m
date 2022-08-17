function h = ineqplot(I,R,c)
% Plotting inequalities, simple and easy.
%
% Input arguments
%            I    -   Inequality as string, i.e.  'x+y>10'
%            R    -   Vector of four components defined by: [xmin, xmax, ymin, ymax], 
%                     if two components are passed: [min, max], the defined region 
%                     will be a square and xmin=ymin=min, xmax=ymax=max.
%            c    -   A three-element RGB vector, or one of the MATLAB 
%                     predefined names, specifying the plot color.
%
% Output arguments
%            h    -   returns the handle of the scattergroup 
%                     object created.
%
% Examples: 
%           >> ineqplot('x.^2+y.^2<10',[-5 5], 'r');
%           >> h = ineqplot('y<x+3',[0 10 -5 5]);
%           >> set(h,'MarkerFaceColor','r'); % Change color
%
% 
% $ Author:   Jorge De Los Santos $
% $ Version:  0.1.0 (initial version) $
%

% Deafault color
if nargin<3, c='b'; end; % blue default
% Length of R vector
if length(R)==4
    xmin = R(1); xmax = R(2);
    ymin = R(3); ymax = R(4);
elseif length(R)==2
    xmin = R(1); xmax = R(2);
    ymin = R(1); ymax = R(2);
end
% hold
set(gca,'NextPlot','add');
% Limits
axis([xmin xmax ymin ymax]);
% Number of points (Change N value for more resolution)
N = 50; % 1000
dx=(xmax-xmin)/N; % 
dy=(ymax-ymin)/N;
[x,y]=meshgrid(xmin:dx:xmax, ymin:dy:ymax);
% Eval condition
idx = find(~eval(I)); 
x(idx) = NaN; %#ok
h = plot(x(:),y(:),'color',c,...
    'LineStyle','-',...
    'MarkerSize',1);
end