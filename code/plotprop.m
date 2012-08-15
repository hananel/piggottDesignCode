function plotprop(filename,Nb,dimension,color,PercentPull,export)
% filename is a csv file with radius, chord and twist columns
% profile is assumed A18 at this point
fignum = 15; figure(15); clf

% read prop properties
M = csvread(filename);
r = M(:,1); c = M(:,2);
Nr = length(r);
c0 = max(c);

% 2D part
if strcmp(dimension,'2D')
    % assuming blade trailing edge runs along radius line
    B.x = [r; r(end:-1:1); r(1)];
    if PercentPull==0           %LE straight
        B.y = [c-c0/2 ;c*0-c0/2 ;c(1)-c0/2];
    else if PercentPull==100    %TE straight - according to first profile
            B.y = [c(1)*ones(Nr,1)-c0/2 ;(c0/2 - c(end:-1:1)) ;c(1)-c0/2];
        else
            % TODO add relation to percentPull for variational pull
            B.y = [c-c0/2 ;c*0-c0/2 ;c(1)-c0/2];
        end
    end
    [B.t,B.r] = cart2pol(B.x,B.y);
    
    % plotting
    figure(fignum); hold on;
    title(['2D properties of' filename]);
    for ii=0:Nb-1
        polar(B.t+ii*(2*pi/Nb),B.r,color);
    end
    
end
% 3D part
if strcmp(dimension,'3D')
    b = M(:,3);
    % read profile x,y
    %Pfilename = '/Users/hananlevy/src/piggott-turbine-design/octave/Input/Section/A18/A18_coordinates.txt';
    [Pfilename, Ppathname, filterindex] = uigetfile({'*.txt';'*.csv'}, 'pick a airfoil coordinate file','Input/Section/GOE417A.txt');
    [P.x0,P.y0,name] = ReadProfile([Ppathname '/' Pfilename]); %Ppathname
    
    figure(fignum+1); clf; hold on;
    for j=0:Nb-1
        for i=1:length(r)
            if PercentPull==0           %LE straight
                [P.t,P.r] = cart2pol(P.x0*c(i),P.y0*c(i));
            else if PercentPull==100    %TE straight - according to first profile
                    if 1
                        % also - correcting TE, to a minimal thickness
                        % specificaly for GOE407A
                        TEmin = 0.6/1000; %m
                        x = P.x0; y = P.y0;
                        xth = x(1:84); yth = y(1:84);
                        xtemp = xth(end); ytemp = yth(end);
                        xtemp(end+1) = 1; ytemp(end+1) = y(1)-TEmin/c(i);
                        
                        % interpolating to fill up reigon with points so that fucking
                        % solidworks doesnt spline it to hell...
                        xth(84:183) = linspace(xtemp(1),xtemp(end),100)';
                        yth(84:183) = interp1(xtemp,ytemp,xth(84:end))';
                        xth(end+1) = 1; yth(end+1) = y(1);
                        % re normilizing
                        P.x0 = xth; P.y0 = yth;
                        [P.t,P.r] = cart2pol(P.x0*c(i)+c0-c(i)-max(c),P.y0*c(i));
                    else
                        [P.t,P.r] = cart2pol(P.x0*c(i)+c0-c(i)-max(c),P.y0*c(i));
                    end
                else
                    % TODO add relation to percentPull for variational pull
                    [P.t,P.r] = cart2pol(P.x0*c(i),P.y0*c(i));
                end
            end
            [P.x,P.y] = pol2cart(P.t+(b(i)*pi/180),P.r);
            [P.t,P.rho,P.zz] = cart2pol(r(i)*ones(length(P.x),1),P.x,P.y);
            [P.x,P.y,P.z] = pol2cart(P.t+j*(2*pi/Nb),P.rho,P.zz);
            if i==1 patch(P.x,P.y,P.z,'k');; end
            plot3(P.x,P.y,P.z,color);
            % saving end points of each blade. trailing edge (1) and
            % leading edge (2)
            Xend1(i) = P.x(1);             Yend1(i) = P.y(1);             Zend1(i) = P.z(1);
            Xend2(i) = P.x(round(end/2));  Yend2(i) = P.y(round(end/2));  Zend2(i) = P.z(round(end/2));
            if export
                csvwrite([filename '_r_' num2str(r(i)) '.txt'],[P.x,P.y,P.z]*1000);
            end
        end
        % plotting contour of blade
        plot3([Xend1 Xend2(end:-1:1)],[Yend1 Yend2(end:-1:1)],[Zend1 Zend2(end:-1:1)],color,'LineWidth',1)
    end
    axis equal
end
