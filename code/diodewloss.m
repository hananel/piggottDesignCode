function  [Pshaft,Pd,RPM]  = diodewloss(inp,plotG);
% 13/8/11
% no loss at the moment
I0 = inp.I0;Kg = inp.Kg; t0 = inp.t0; Rw = inp.Rw; v0 = inp.v0; a = inp.RPMmin; b = inp.RPMmax;


% Arranging generator parameters (specifically for DC generator as used in Adital Ela's design)
Rtot = Rw; %[ohm]
RPM = a:10:b;

V =RPM*Kg;
I = diode_IV(V,I0,t0,v0); 

Pd = (V-v0).* I;
Pshaft = I.*Kg.*RPM; % old model was Pshaft = Vd*I + I.^2*Rtot; 

eff = Pd./Pshaft*100;

if nargin == 2
    f = figure(10); hold on;
    set(f,'Position',[000 000 500 300 ])
    Title('Generator model');
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
    legend([H3,H1,H2],'Power out', 'Power in', 'Efficiency','Location','SouthEast')
    set(AX(2),'Ylim',[0 100]); set(AX(2),'YTick',0:10:100);
    ylabel(AX(1),'Power [W]'); ylabel(AX(2),'Efficiency [%]'); xlabel('RPM');
    
    % making png and svg
%    plot2svg([inp.CpCt_fileDirectory '/Gen_' inp.filename '.svg']);
    print( gcf, '-dpng', ['../' inp.CpCt_fileDirectory '/Gen_' inp.filename '.png'])
end
