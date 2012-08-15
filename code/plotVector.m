function plotVector(XYs,XYe,col,options)
% XYs is the start coordinates of the vector
% XYe is the end coordinates of the vector
hold on;
if nargin<4
    handles = plot_arrow( XYs(1) ,XYs(2), XYe(1) ,XYe(2) ,'color',col,'facecolor',col,'edgecolor',col);
else
    handles = plot_arrow( XYs(1) ,XYs(2), XYe(1) ,XYe(2) ,'color',col,'facecolor',col,'edgecolor',col,options);
end

