function [minimum] = PowerCurve(inp,constantPMG)
% Hanan E. Levy 4/2010
% minimum - for use with
% fminbnd(@(x) PowerCurve(x),4,40)
% gives 10^5 W minus the total power summation under the power curve (sum(P))

%load InpFile

% Input_NewBlade.m reads inp.file_blade geometry and interpolates to inp.Nr
% elements

inp.normalRun = 1;
Input_NewBlade

WorldConstants

Vvec = inp.Vfstart:inp.VvecStep:inp.Vfend;
firstTime = 1;

if(~inp.plotoff)
    f = figure(1); subplot(321);
    set(f,'Position',[500 500 500 400 ])
end
i = 0;

% input - from seperate file
[CpData CtData TSR VData] = loadCpCtData(inp);


maxCp = max(max(CpData(2:end,2:end)));
[row,col]=find(CpData==maxCp,1);
TSRopt = TSR(row);
if ~inp.measuredPMG
    [inp.Kg,inp.Rg,inp.coil] = PMGdesign(inp.Rstator,inp.magnet,inp.coil);
end

if (~inp.plotoff && inp.plotGen)
    [inp.coil, inp.magnet] = plotStator(inp.coil,inp.magnet,inp.Rstator*1000,inp.Rw,inp.r,inp.measuredPMG);
    set(gcf,'Color','w');
    %plot2svg([inp.file_Output_directory '/startor_' inp.filename '.svg']);
    print( gcf, '-dpng', [inp.file_Output_directory '/startor_' inp.filename '.png'])
end

inp.RPMmax = TSR(end) * Vvec(end) / R * 30/pi*1.1;         %[RPM] max range of generator
inp.RPMmin = 0;
i = 0;

if ~inp.plotoff figure(1); hold on;end

for V = Vvec
    i = i+1;
    % battery changing voltage with extra current - data from measurement
    % of CometME 4.2 meter 48 v machine
    if inp.loadVlaw
        load BattLaw_CometME42
        inp.Vb = interp1(BattLaw(:,1),BattLaw(:,2),V);
    end
    % turbine power curve
    RPMt = TSR * V / R * 30/pi;     %[RPM]
    if length(VData)==1
        Cp = interp1(TSR,CpData(2:end,2:end),TSR)*rho/rho0;
    else
        Cp = interp2(VData,TSR,CpData(2:end,2:end),V,TSR)*rho/rho0;
    end
    Pt = 0.5*rho*V^3*pi*R^2.*Cp;    %[W]
    % generator power curves
    % if using voltge control
    if inp.useVoltageControl
        inp.Vb = interp1(inp.controlSpeed,inp.controlVolt,V)*(inp.VoltageControlDelta+1);
    end
 
    if ~inp.plotoff && inp.plotGen && firstTime
        [Pg,Pbatt,RPMg] = PMGwloss(inp,1);
        firstTime = 0;
    else
        if inp.diode
            [Pg,Pbatt,RPMg] = diodewloss(inp);
        else
            [Pg,Pbatt,RPMg] = PMGwloss(inp);
        end
    end
    
    if ~inp.plotoff figure(1); subplot(321); hold on;end
    Peff = Pbatt./Pg;
    for j=1:length(Peff)
        if isnan(Peff(j))
            Peff(j) = 1;
        end
    end
    % Constant speed generator - kept for demonstration purposes
    if nargin==2  % constant speed generator
        RPMg = [constantPMG constantPMG]; Pg = [0 max(Pt)*1.1];
        Rtemp = constantPMG;  %[RPM]
        Ptemp = interp1(RPMt,Pt,constantPMG);
    else                    % PMG (with loses)
        [Rtemp,Ptemp]=polyxpoly(RPMg,Pg,RPMt,Pt);
    end
    if ~isempty(Rtemp)
        RPM(i) = Rtemp(1);
        P(i) = Ptemp(1);
        if length(VData)==1
            Ct = interp1(TSR,CtData(2:end,2:end),TSR);
        else
            Ct = interp2(VData,TSR,CtData(2:end,2:end),V,TSR);
        end
        T(i) = 0.5*rho*V^2*pi*R^2.*interp1(TSR,Ct,RPM(i)*pi/30*R/V);
        
        %ploting
        if(~inp.plotoff)
            if (nargin==2)
                plot(RPMt,Pt,'b',RPMg,Pg,'k',RPM(i),P(i),'*r');
            else
                Geff = interp1(RPMg,Peff,RPM(i));
                Pcharge(i) = Geff*P(i);
                plot(RPMt,Pt,'b',RPMg,Pg,'k',RPMg,Pbatt,'k:',RPM(i),P(i),'*r');
            end
        else
            Geff = interp1(RPMg,Peff,RPM(i));
            Pcharge(i) = Geff*P(i);
        end
        if ~inp.plotoff axis([0 max(RPMt) 0 max(Pt)]); end
        
        % calculating local reynolds at 0.75 radius
        Vbe = sqrt(V^2 + (0.75*R*pi/30*RPM(i))^2);
        Re(i) = rho/ni*Vbe*inp.chord75;
        VRe(i) = V;
    end
end

% if no match is found
if ~exist('RPM')
    disp('no solution with this generator')
    minimum = 0;
    return
end
% if turbine stalls at high wind speed
if length(RPM)<length(Vvec)
    RPM(end+1:length(Vvec)) = 0;
    P(end+1:length(Vvec)) = 0;
    Pcharge(end+1:length(Vvec)) = 0;
end
% calculating TSR
TSRi = RPM .* pi/30 * R ./ Vvec;
% calculating efficiency
Beff = P ./ (0.5*rho*Vvec.^3*pi*R^2)*100; % [%]
if ~(nargin==2)
    Geff = Pcharge ./ ( 0.5*rho*Vvec.^3*pi*R^2)*100; %[precent]
end
if(~inp.plotoff)
    figure(1); subplot(321); hold on;
    xlabel('RPM'); ylabel('P [Watt]')
    if ~(nargin==2)
        
        if inp.compare2TestData         % comparing to data
            TestData = csvread([inp.TestData_fileDirectory '/' inp.TestData_file]);
            %TestData(:,1) = TestData(:,1).*(12.5/10)^(0.1);
            % efficiency plots
            if mean(TestData(:,2))<1
                subplot(3,2,2);hold on;
                plot(TestData(:,1),0.5*rho*pi*inp.R^2.*TestData(:,1).^3.*TestData(:,2),'r.') %
                subplot(3,2,[3 4]);hold on;
                plot(TestData(:,1),100*TestData(:,2),'r.');
            else
                subplot(3,2,2); hold on;
                plot(TestData(:,1).*(12.5/10)^(0.14),TestData(:,2),'r.')
                subplot(3,2,[3 4]); hold on;
                plot(TestData(:,1),100*TestData(:,2)./(.5*rho*TestData(:,1).^3*pi*inp.R^2),'r.');
            end
            axis([Vvec(1) Vvec(end) 0 59.3])
            plot(Vvec,Beff,inp.color); plot(Vvec,Geff,inp.color,'LineWidth',2);
            legend('Data [%]','Blade eff [%]','Overall eff [%]');
            ylabel('Efficiency [%]')
            subplot(3,2,[5 6]); hold on;
            plot(Vvec,TSRi,inp.color,'LineWidth',2);
            plot([min(Vvec) max(Vvec)],[TSRopt TSRopt],inp.color);
            if size(TestData,2)>2
                plot(TestData(:,1),TestData(:,3),'r.');
                legend('TSR','TSRopt','Data');
            else
                legend('TSR','TSRopt');
            end
        else
            subplot(3,2,[3 4]);hold on;
            hold on; plot(Vvec,Beff,inp.color); plot(Vvec,Geff,inp.color,'LineWidth',2);
            legend('Blade eff [%]','Overall eff[%]');
            ylabel('Efficiency [%]')
            subplot(3,2,[5 6]); hold on;
            plot(Vvec,TSRi,[ inp.color ],[min(Vvec) max(Vvec)],[TSRopt TSRopt],inp.color);
            legend('TSR','TSRopt');
        end
        
        xlabel('V [m/s]'); ylabel('TSR')
    end
    subplot(322); hold on;
    if (nargin==2)
        plot(Vvec,P/1000);
    else
        plot(Vvec,Pcharge,inp.color);
        ylabel('P in battery [W]')
        xlabel('V [m/s]');
        grid on;
    end
end
V = Vvec;

% Reynolds range
% figure(20); hold on; plot(VRe,Re); title('Reynolds at 0.75R');

% TODO: missing - thrust range!

% output
if (~inp.plotoff)
    figure(1);
    set(gcf,'Color','w');
end
% plot2svg([inp.file_Output_directory '/PowerCurve_' inp.filename '.svg']);
% print( gcf, '-dpng', [inp.file_Output_directory '/PowerCurve_' inp.filename '.png'])
if length(Vvec)>1
    % optimizer
    % maximum energy
    WBLscale = 5; WBLshape = 2;
    windPdfInRange = exp(-(V(1:end-1)./WBLscale).^WBLshape) - exp(-(V(2:end)./WBLscale).^WBLshape);
    dV = V(2)-Vvec(1); Venergy = V(1:end-1)+dV/2;
    energyDist = 0.5*(Pcharge(1:end-1) + Pcharge(2:end)).*windPdfInRange * 24; %[kwh/day / (dV) ]
    minimum = -abs(sum(energyDist));
    disp(['For a hub height condition of c,k = ' num2str(WBLscale) ',' num2str(WBLshape) ', output is: ' num2str(sum(energyDist)) ' wh'])
    
    %available energy
    P = 0.5*rho*Vvec.^3*pi*R^2;
    Edist = 0.5*(P(1:end-1) + P(2:end)).*windPdfInRange * 24; %[wh/day / (dV) ]
    E = sum(Edist);
    eff = -minimum/E;
    
    %plotting energy dist
    if(~inp.plotoff)
        f = figure(3); hold on;
        set(f,'Position',[500 000 500 300 ])
        [AX,H1,H2] = plotyy(Venergy,energyDist,Venergy,windPdfInRange*24);
        set(get(AX(1),'Ylabel'),'String','Wh/day / (m/s)')
        set(get(AX(2),'Ylabel'),'String','hours per day')
        xlabel('V m/s');
    end
else
    disp('WindyBoy optimization point:');
    disp(['V[m/s] = ' num2str(Vvec) ' Pcharge[Watt] = ' num2str(Pcharge) ' Volts = ' num2str(inp.Vb)]);
    minimum = 0;
end
%fflush(stdout);

%load SWTpowercurves
%plot(BWC850(:,1),BWC850(:,2),'r',marlec(:,1),marlec(:,2),'b',Air403(:,1),Air403(:,2),'g',AA100(:,1),AA100(:,2),'y')

