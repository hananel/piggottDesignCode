function   [minimum] = OptimizeNewBlade(parameter, inp)
% Inputs:
%  inp.TSRstart - starting tip speed ratio
%  inp.TSRend   - ending tip speed ratio
%  inp.beta0    - pitch angle correction term
%  inp.doplots     - 1 for displaying graphs, 0 for not
% use for optimization: fminbnd(@(parameter) OptimizeNewBlade(parameter,inp),0 ,5)
%
% OptimizeNewBlade uses BEM therory to calculate the Cp and Ct vs. TSR of a blade
% geometry, giving in the program by chord and Beta, and Cl/Cd as a
% function of alpha are given in ActuatorDisk.m
%
%     Cmax          = The maximum chord length
%     Cl ,Cd        = 2-D lift and drag coefficients, respectively
%     Ct ,Cp        = Thrust and power coefficients, respectively
%     e             = Cross-sectional drag to lift ratio
%     Nb            = Number of blades
%     Pa            = Ambient pressure
%     dP            = Power of a disk annulus
%     r             = Radial coordinate at the disk plane
%     R             = Turbine radius
%     Re            = Reynolds number
%     dT            = Thrust increment of disk annulus
%     Vbe           = Resultant cross-sectional velocity
%     Vf            = Wind speed
%     aa            = Axial induced velocity ratio components, at the disk plane
%     at            = circumferential induced velocity ratio components, at the disk plane
%     alpha         = Angle of attack
%     beta          = Pitch angle distribution along the blade
%     phi           = Inflow angle
%     TSR           = Tip speed ratio
%     tsr           = Local tangential speed ratio
%     W             = Turbines rotational speed

% declarations

beta = 0;
inp_Windylight_normal_fucking_turbine
filename = inp.filename;
doplots = inp.doplots;


% Input_NewBlade.m reads inp.file_blade geometry and interpolates to inp.Nr
% elements
Input_NewBlade
% use parameter "optimizeParameter" as parameter for optimization
% TODO insert code that rewrites a parameter acording to the input in optimizeParameter
% something like
% switch inp.optimizaParamName
%   case 'c'
%       blade.c = optimizaeParameter
%     .
%     .
%     .
%

% this is the shape taken - the optimal distibution.
% The parameter being varied is TSR, andthe claulation looks for maximum
% energy over typical wind speeds
blade.Nb = 6;
% for optimization criterion -sum(max(Cp')); for V = 3:8, TSR = 0.5:0.5:5
% Nb = 4, TSR_opt = 3.6152 
% Nb = 5, TSR_opt = 3.5497
% Nb = 6, TSR_opt = 3.6449 

[blade.r,blade.c,blade.beta] = optimal_blade_no_drag(parameter,R,0.2*R,blade.Nb,inp.Nr,0.055,1.0,8*pi/180,'A18','k',1,1);
blade.dr = blade.r(2)-blade.r(1);

disp([ ' blade number = ' num2str(blade.Nb)])
j=0;
% loop through TSR range
VfVec = inp.Vfstart:inp.Vfstep:inp.Vfend;
% color coding
Vl = length(VfVec);
color = jet(Vl);
%beta = blade.beta*180/pi;
for Vf = VfVec
    i=0;
    
    j = j+1; inp.Vf = Vf;
    TSRvec = inp.TSRstart:inp.TSRstep:inp.TSRend;
    fprintf(1,'\nVf = %2.1f m/s',Vf);
    for TSR = TSRvec
        fprintf(1,'\n   TSR = %4.2f ',TSR); %fflush(stdout);
        inp.TSR = TSR;
        i = i+1;
        if inp.useOhad % Ohad Gur 2008 paper implemented by Mickey Reitman
            BE_M_method(i) = ActuatorDisc(inp,blade);
            dCt   = real(killNans(BE_M_method(i).dT./(rho*Vf^2*pi.*r*dr),1));
            dCp   = real(killNans(BE_M_method(i).dP./(rho*Vf.^3*pi.*r*dr),1));
            %Ct(j,i) = sum(dCt)*dr/R;
            P(i)  = sum(BE_M_method(i).dP);
            %Cp(j,i) = sum(dCp)*dr/R;              % TODO verify
            T(i)  = sum(BE_M_method(i).dT);
            M(i)  = sum(BE_M_method(i).dT.*(r-r(1)));
            Cp(j,i) = P(i)/(0.5*rho*Vf^3*pi*(r(end)+blade.dr)^2);
            Ct(j,i) = T(i)/(0.5*rho*Vf^2*pi*(r(end)+blade.dr)^2);
        else
            %Danish book "wind turbine aerodynamics" 2000 implemented by Hanan einav Levy
            BE_M_method(i) = ActuatorDiscH(inp,blade);
            dCt   = BE_M_method(i).dT./(rho*Vf^2*pi.*r*dr);
            dCp   = BE_M_method(i).dP./(rho*Vf.^3*pi.*r*dr);
            %Ct(j,i) = sum(dCt)*dr/R;
            %Cp(j,i) = sum(dCp)*dr/R;              % TODO verify
            P(i) = sum(BE_M_method(i).dP);
            T(i) = sum(BE_M_method(i).dT);
            M(i)  = sum(BE_M_method(i).dT.*(r-r(1)));
            Cp(j,i) = P(i)/(0.5*rho*Vf^3*pi*(r(end)+blade.dr)^2);
            Ct(j,i) = T(i)/(0.5*rho*Vf^2*pi*(r(end)+blade.dr)^2);
        end
        
        inp.aa = BE_M_method(i).aa;
        inp.at = BE_M_method(i).at;
        
    end
    if (inp.doplotsCp)
        colorPos = j;
        figure(2); hold on;
        subplot(211);
        hold on;  
        if inp.useOhad
            hCp = plot(TSRvec,Cp(j,:),'-','Color',color(colorPos,:));
            grid on;
        else
            hCp = plot(TSRvec,Cp(j,:),'-.','Color',color(colorPos,:));
            grid on;
        end
        title('Power and force coeffieicnts','FontSize',16); ylabel('Cp','FontSize',16);  axis tight
        subplot(212);
        hold on;  
        if inp.useOhad
            hCt = plot(TSRvec,Ct(j,:),'-','Color',color(colorPos,:));
            grid on;
        else
            hCt = plot(TSRvec,Ct(j,:),'-.','Color',color(colorPos,:));
            grid on;
        end
        ylabel('Ct','FontSize',16);axis tight
        xlabel('TSR','FontSize',16);
        set(gcf,'Color','w')
    end
    
end

% hard output
csvwrite(['../' inp.CpCt_fileDirectory inp.Cp_filename], [0 VfVec ;TSRvec' real(Cp)'])
csvwrite(['../' inp.CpCt_fileDirectory inp.Ct_filename], [0 VfVec ;TSRvec' real(Ct)'])
print( gcf, ['../' inp.CpCt_fileDirectory '/' inp.CpCt_filename '.png'], '-dpng')
fprintf(1,'\n');
disp(['wrote output to ../' inp.CpCt_fileDirectory '/' inp.Cp_filename])
disp(['            and ../' inp.CpCt_fileDirectory '/' inp.Ct_filename])

% minimization criterion
minimum = -sum(max(Cp'));

% TODO
% insert a change to the inp.file_CpCt to the current location of the new
% output file!
% matlab: plot2svg([inp.CpCt_fileDirectory '/' filename '.svg'],gcf,'png');
% matlab: print( gcf, '-dpng', [inp.CpCt_fileDirectory '/' inp.CpCt_filename '.png'])
%g = get(gcf); g = g.children; set(g(1),'FontSize',14); set(g(2),'FontSize',14);


