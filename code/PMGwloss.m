function  [Pshaft,Pb,RPM]  = PMGwloss(inp,plotG); % eff_500

Vb = inp.Vb;Kg = inp.Kg; Rg = inp.Rg; Rw = inp.Rw; r = inp.r; a = inp.RPMmin; b = inp.RPMmax;

% a = RPM
% a = RPMfirst, b=RPMlast
% generator charestaristics (from experiment)
%RPMcut    %[RPM]       when battery voltage is accesed
%RPMslope  %[W / RPM]   ratio of watts per RPM increase, after cut-in
%
% from the pico turbine pdf:
% V = N*A*R*P*B/2 = Kg * RPM = Kg * R * 60 ==> Kg = N*A*P*B/120
% N = number of loops of wire WindAid 26
% A is the area enclosed by a loop of wire, in square meters. 40/100^2
% R is the rotational velocity of the magnets, in cycles per second.
% R = RPM/60
% P is the number of magnet poles per cycle. WindAid: 24
% B is the strength of the magnetic field of each pole, in Tesla. ~1
% P = I*Vb = (kg*RPM - Vb)/Rr*Vb
%
%A = 0.005; % [m^2]
%B = 1.2;     % [T]
%N = 26;
%P = 48;
%R = RPM/60;

% Arranging generator parameters
Vd = 1.4;                 %[V] diode voltage drop
Rtot = Rw + Rg; %[ohm]
RPMcut = (Vb+Vd)/Kg; % Vb = 24; Kg = 24/100
Powerslope = Kg*Vb/Rtot;% Rr = 0.288

if a>RPMcut
    RPM = a:10:b;
else
    RPM = RPMcut:10:b; %[a RPMcut:10:b]
end

if RPMcut>b
    disp('RPMcut > b')
    return
end

I = (RPM*Kg-Vb-Vd)./(Rtot+r); %[A]
condition = 1;
A = 0; %0.0001;
R0 = Rtot;
I0 = I;
while condition
    % correcting for increase in resistence due to current -
    % taking Rtot to be Rg. for simplification
    Inew = I;
    Rtot = R0.*(1+I.^2*A);
    Itag = (RPM*Kg-Vb-Vd)./(Rtot+r); %[A]
    I = Itag - 0.3*(Itag-I); % relaxation
    convergence = abs(sum((Inew-I)./I));
    if convergence<0.01, condition=0; end
end

if inp.WindyBoy
    V = Vb;
else
    V =RPM*Kg -I.*Rtot;
    for i=1:length(V)
        if V(i)>31/24*Vb
            V(i) = 31/24*Vb;
        end
    end
end
Pb = V.* I;
Pshaft = I.*Kg.*RPM; % old model was Pshaft = Vd*I + I.^2*Rtot;
for i=1:length(Pb)
    if (Pb(i)<0) || (I(i)<0)  Pb(i)=0;
    end
    if (Pshaft(i)<0) || (I(i)<0) Pshaft(i)=0;
    end
end
eff = Pb./Pshaft*100;
eff_1200 = eff(find(Pb>1200,1));

% empirical correction - 
RPM = RPM .* linspace(1,inp.RPMcor,length(RPM));
if nargin == 2
    f = figure(10); hold on;
    set(f,'Position',[000 000 500 300 ])
    title({'Generator model','Efficiency in yellow','Shaft power in magenta','Battery power in blue'});
    RPMminPos = find(Pshaft>0,1,'first');
    RPMmaxPos = find(eff>40,1,'last');
    RPMplot = RPM(RPMminPos:RPMmaxPos);
    [AX,H1,H2] = plotyy(RPMplot,Pshaft(RPMminPos:RPMmaxPos),RPMplot,eff(RPMminPos:RPMmaxPos)); hold on;
    set(H1,'Color','k');            set(H2,'Color','k')
    set(H1,'Marker','s');           set(H2,'Marker','^');
    set(H1,'MarkerFaceColor','m');  set(H2,'MarkerFaceColor','y');
    H3 = plot(AX(1),RPMplot,Pb(RPMminPos:RPMmaxPos),'kd-','MarkerFaceColor','b');
    set(AX(1),'XGrid','on')
    set(AX(1),'YGrid','on')
    % plotting voltage as text along power out axis
    RPMtextNum = 8;
    RPMspacedVec = round(linspace(RPMminPos,RPMmaxPos,RPMtextNum));
    shiftx = 2*(RPMspacedVec(2)-RPMspacedVec(1));
    shifty = (Pb(end)-Pb(1))/(RPMtextNum*10);
    for i=RPMspacedVec(2:end-1)
        text(RPM(i)+shiftx,Pb(i)-shifty,[num2str(V(i),'%3.1f') ' V']);
    end
    %legend([H3,H1,H2],'Power out', 'Power in', 'Efficiency','Location','SouthEast')
    set(AX(2),'Ylim',[0 100]); set(AX(2),'YTick',0:10:100);
    ylabel(AX(1),'Power [W]'); ylabel(AX(2),'Efficiency [%]'); xlabel('RPM');
    
    % making png and svg
    %    plot2svg([inp.CpCt_fileDirectory '/Gen_' inp.filename '.svg']);
    print( gcf, '-dpng', ['../' inp.CpCt_fileDirectory '/Gen_' inp.filename '.png'])
end
