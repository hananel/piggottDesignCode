function [coil, magnet] = plotStator(coil,magnet,R,Rw,r,measuredPMG)

% position of plot
g = figure(2);
set(g,'Position',[000 500 500 400 ])

% numerical parameters
N = 50;
t = 10;
f = 1.3;

% shape parameters
x = linspace(0,coil.legWidth,N);
y = sqrt(coil.legWidth^2-x.^2);
xr = linspace(coil.legWidth,0,N);
yr = sqrt(coil.legWidth^2-xr.^2);

% counturs
coil.x = [magnet.width/2+x,magnet.width/2+coil.legWidth-x,-magnet.width/2-x,-magnet.width/2-xr,magnet.width/2];
coil.y = [R+magnet.length/2+y,R-magnet.length/2-yr,R-magnet.length/2-y,R+magnet.length/2+yr,R+magnet.length/2+coil.legWidth];
magnet.x = [magnet.width/2,magnet.width/2,-magnet.width/2,-magnet.width/2,magnet.width/2];
magnet.y = [R+magnet.length/2,R-magnet.length/2,R-magnet.length/2,R+magnet.length/2,R+magnet.length/2];
[coil.t,coil.r] = cart2pol(coil.x,coil.y);
[magnet.t,magnet.r] = cart2pol(magnet.x,magnet.y);

magnet.alpha = 2*pi/magnet.N;
coil.alpha = 2*pi/coil.N;
 
subplot(1,2,1); hold on;
for i=0:coil.N-1
    polar(coil.t+i*coil.alpha,coil.r,'r');
end
for i=0:magnet.N-1
    polar(magnet.t+i*magnet.alpha,magnet.r,'b');
end

% added - 29/10/12 Hanan Einav-Levy
% Drawing B and C limiting coil circles - according to Hugh Piggott plans (2012) page 40 
C = floor(min(coil.r));
B = ceil(max(coil.r));
circle = 0:0.01:2*pi;
polar(circle,B*ones(1,length(circle)),'k')
polar(circle,C*ones(1,length(circle)),'k')
g = line([0,0],[0 R]);
set(g,'LineStyle','-.','Color','k','Marker','x');
text(t,R/2,num2str(R));
axis([-(R+magnet.length/2)*f (R+magnet.length/2)*f -(R+magnet.length/2)*f (R+magnet.length/2)*f])
axis off
axis square
xlabel('mm')
disp('')
if measuredPMG
    titleDesignProperties = sprintf('Design properties\nused just for drawing\nPMG data from measurements');
    disp(titleDesignProperties)
else
    if coil.AWG
        titleDesignProperties = sprintf('Design properties\nMagnet field strength = %2.2f Tesla\nFF = %2.0f%%\nAWG %2.0f# X %1d with %2.0f turns',magnet.B,coil.FF*100,coil.AWG,coil.wireN,coil.turns);
    else
        titleDesignProperties = sprintf('Design properties\nMagnet field strength = %2.2f Tesla\nFF = %2.0f%%\nCoil %2.1f mm diameter X %1d with %2.0f turns',magnet.B,coil.FF*100,coil.mm, coil.wireN,coil.turns);
    end
    textStr = {titleDesignProperties,['p. 40 limiting values of B = ',num2str(B),' [mm], C = ',num2str(C),' [mm]'],['Coil thickness (maximum mold thickness) = ',num2str(coil.thickness),' [mm]'],['Coil number = ',num2str(coil.N)],['Magnet number = ',num2str(magnet.N)]};
    for ii=1:length(textStr)
        disp(textStr{ii})
    end
    text(-B,2*B,textStr)
end

subplot(1,2,2); hold on;
plot(coil.x,coil.y-R,'r');
plot(magnet.x,magnet.y-R,'b');
g = line([magnet.width/2,magnet.width/2+coil.legWidth],[0 0]);
set(g,'LineStyle','-.','Color','k','Marker','x');
text((magnet.width/2+coil.legWidth)/2,t/3,num2str(coil.legWidth,3));
text(-magnet.width,0,num2str(magnet.length,3));
text(-t/2,magnet.length/2+t/2,num2str(magnet.width,3));
axis equal
axis off
axis square
xlabel('mm')

if measuredPMG
    titleMeasuredPMG=sprintf('PMG properties\nKg = %2.1f [volt/RPM]\nRg =  %2.1f\n\nSystem properties\nRcable = %2.2f Ohm\nRbatt = %2.2f Ohm',coil.Kg,coil.Rg,Rw,r);
else
    titleMeasuredPMG=sprintf('PMG properties\nKg = %2.1f [volt/RPM]\nRg = %2.1f Ohm-impedence\n\nSystem properties\nRcable = %2.2f Ohm\nRbatt = %2.2f Ohm',coil.Kg,coil.Rg,Rw,r);
end
disp(titleMeasuredPMG)
subplot(1,2,1); text(-B,-B*1.6,titleMeasuredPMG)
