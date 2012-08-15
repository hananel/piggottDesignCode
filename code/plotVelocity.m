function plotVelocity(Vf,TSR,r,offset,col)
% Vf = wind velocity [m/s]
% TSR = Tip Speed Ratio
% r = ratio of plot size to Vf
% offset = offset of center of forces (positive in the negative x direction)
% col = color of velocity vectors

if exist('octave_core_file_name')
    %doesn't work on octave
    return
end

hold on;
x = 8;
plot_arrow( -offset ,-Vf/r, -offset,0 ,'color',col,'facecolor',col,'edgecolor',col);
plot_arrow( -(Vf*TSR/r+offset) ,0, -offset,0 ,'color',col,'facecolor',col,'edgecolor',col,'headwidth',0.07/x,'headheight',0.15/x);
plot_arrow( -(Vf*TSR/r +offset), -Vf/r, -offset,0 ,'facecolor','k','edgecolor','k','headwidth',0.07/x,'headheight',0.15/x*4);