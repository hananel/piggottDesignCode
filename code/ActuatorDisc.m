function  BE_M_method =  ActuatorDisc(inp,blade)
% BEM actuator disk model.
% complementry to NewBlade.m
% uses NACA4412.m to calculate Cl and Cd as a function of alpha, according to table
% no reynolds related changes where accounted for at this time
% 30.8.10

% MultiStart('UseParallel', 'always');

Vf = inp.Vf;
TSR = inp.TSR;
draw = inp.doplots;
iter = inp.iter;
maxiter = inp.maxiter;
relax = inp.relax;

Nb = blade.Nb;
R = blade.Dm/2;
r0 = blade.r(1);
r = blade.r;
dr = blade.dr;
beta = blade.beta;
chord = blade.c;

%WORLD CONSTANTS
WorldConstants

Nr = length(r);              % number of blade elements
vtip = Vf*TSR;
W  = vtip/(R);               % [radv/sec] rotational speed

if isfield(inp,'aa')
    wa = inp.aa;
    wt = inp.at;
else
    wa = 0.05*ones(Nr,1);
    wt = 0.05*r.*ones(Nr,1);
end

convergence_criteria1 = 100;
convergence_criteria2 = 100;
iter_criteria = 1;
epsilon = inp.epsilon;
wa_old = 100; wt_old = 100;
if TSR>7 relax = relax/2; end
while and(or(max(convergence_criteria1)>epsilon,max(convergence_criteria2)>epsilon),iter_criteria)
    
    fprintf(1,'.'); %fflush(stdout);
    %sinPhi = (Vf-wa)./Vbe;
    Vbe = sqrt( (Vf-wa).^2 + (W*r+wt).^2);
    sinPhi = min((Vf-wa)./Vbe,0.95);
    phi = asin(sinPhi);
    cosPhi = cos(phi);
    alpha = phi - beta;
    alphaDeg = alpha / pi*180;
    Re = rho/ni.*Vbe.*chord;
    
    % TODO: aeroelasticc correction
    % what we need:
    %
    
    % calculating blade element properties
    [Cl Cd] = LoadSectionS(real(alphaDeg), real(Re), blade); 
    
    Cd0 = inp.Cd0; % this is the fudge factor for adital's high solidity thick unaerodynamic design
    
    dT = 0.5*rho.*Vbe.^2*Nb.*chord.*(Cl.*cosPhi+Cd.*sinPhi)*dr;             %eq. (7)
    dP = 0.5*rho.*Vbe.^2*Nb.*chord.*(Cl.*sinPhi-Cd.*cosPhi)*W.*r*dr  ;        %eq. (8)
    dT = killNans(dT);  dP = killNans(dP);
    % MICKY - where is this from ???
    dTlimit = 4*pi*rho*r*Vf^2*dr/4;
    dT = min(dT, dTlimit);
    
    k = Gold(Nb*ones(Nr, 1), sinPhi, r./R);k = real(k);
    aa = real((1-sqrt(max(1-(dT./dr)./(pi*rho.*r*Vf^2),0)))./(2.*k));       %eq. (2) solution
    at = 1./(4*pi*rho.*r.^3*Vf*W^2.*k.*(1-k.*aa)).*(dP./dr);                %eq. (3) solution
    delta_a = wa - aa*Vf;
    delta_t = wt - at*W.*r;
    
    % relaxation
    wa = wa - relax*delta_a;
    wt = wt - relax*delta_t;
    
    %convergence
    convergence_criteria1 = abs((wa-wa_old)./wa);
    convergence_criteria2 = abs((wt-wt_old)./wt);
    
    wt_old = wt;
    wa_old = wa;
    if iter>maxiter
        iter_criteria=0;
    end
    iter = iter+1;
    
    if draw
        f = figure(1); set(f,'Position',[50 100 600 500 ])
        subplot(3,2,3);
        title('red: w_t , black: w\_a ');
        plot(r,wt,'r',r,wa,'k');hold off
        axis tight;
        title('induced wind speed');
        
        dCp = dP./(rho*Vf.^3*pi.*r*dr);
        dCt = dT./(rho*Vf^2*pi.*r*dr);
        subplot(3,2,[1 2]);
        plot(r,dP,'-r',r,dT,'k');legend('dP','dT','Location','northeast');
        %axis([0, r(end), min(dCp) ,max(dCp)]);
        title(['Forces along blade TSR = ' num2str(TSR) ' V = ' num2str(Vf) ' m/s iteration #' num2str(iter)]);
        grid on;
        
        subplot(3,2,[4 6]);
        plot(r,alphaDeg,'.');
        title('Blade angles');   hold on
        plot(r,beta*180/pi,'-r', r,phi*180/pi,'g');hold off
        %axis([0 max(r) min([alphaDeg beta*180*pi]) max([alphaDeg beta*180*pi])]);
        grid on;
        
        subplot(3,2,5);
        plotyy(r,Cl./Cd,r,Re);
        title('Cl/Cd and Re');
        grid on;
        
        % experimental - plotting velocity and forces at 4 cross
        % TODO - ADD FORCE PLOTTING
        % sections along blade radius
        if 0==1
            gg = figure(10); clf; set(gg,'Position',[670 100 500 500 ]); sc = 1.2; %scale parameter for axis size
            % choosing cross sections
            r1 = 1;                                                 %root
            subplot(221);
            plotNACAXXXX([4400 + thick(r1)],chord(r1),180/pi*beta(r1),r(r1),'k');
            plotVelocity(Vf,TSR*r(r1)/R,10/chord(r1),0.25*chord(r1),'r') ; axis([-chord(1)*sc chord(1)*sc -chord(1)*sc chord(1)*sc]);
            r2 = round(length(r)/3);            % 1/3 radius
            subplot(222);
            plotNACAXXXX([4400 + thick(r2)],chord(r2),180/pi*beta(r2),r(r2),'k');
            plotVelocity(Vf,TSR*r(r2)/R,10/chord(r1),0.25*chord(r2),'r'); axis([-chord(1)*sc chord(1)*sc -chord(1)*sc chord(1)*sc]);
            r3 = round(length(r)/3*2);         % 2/3 radius
            subplot(223);
            plotNACAXXXX([4400 + thick(r3)],chord(r3),180/pi*beta(r3),r(r3),'k');
            plotVelocity(Vf,TSR*r(r3)/R,10/chord(r1),0.25*chord(r3),'r'); axis([-chord(1)*sc chord(1)*sc -chord(1)*sc chord(1)*sc]);
            r4 = length(r);                                % tip
            subplot(224);
            plotNACAXXXX([4400 + thick(r4)],chord(r4),180/pi*beta(r1),r(r4),'k');
            plotVelocity(Vf,TSR*r(r4)/R,10/chord(r1),0.25*chord(r4),'r'); axis([-chord(1)*sc chord(1)*sc -chord(1)*sc chord(1)*sc]);
            %drawnow("expose")
        end
        drawnow
    end
    % detecting non convergence hysteresis due to bad 2D data
    % so simply checking if every 2nd iteration is the same, under epsilon
    % difference
    %saving results vector
    dPmat(iter,:) = dP;
    dTmat(iter,:) = dT;
    phimat(iter,:) = phi;
    alphamat(iter,:) = alpha;
    Remat(iter,:) = Re;
    aamat(iter,:) = aa;
    atmat(iter,:) = at;
    if iter>3
        if (abs(sum(dPmat(end,:) - dPmat(end-2,:))./mean(dPmat(end,:)))<epsilon)
            iter_criteria=0;
            dP = (dPmat(iter,:)' + dPmat(iter-1,:)')/2;
            dT = (dTmat(iter,:)' + dTmat(iter-1,:)')/2;
            phi = (phimat(iter,:)' + phimat(iter-1,:)')/2;
            alpha = (alphamat(iter,:)' + alphamat(iter-1,:)')/2;
            Re = (Remat(iter,:)' + Remat(iter-1,:)')/2;
            aa = (aamat(iter,:)' + aamat(iter-1,:)')/2;
            at = (atmat(iter,:)' + atmat(iter-1,:)')/2;
        end
    end
end

dT = killNans(dT);  dP = killNans(dP);
BE_M_method.dP = dP;
BE_M_method.dT = dT;
BE_M_method.phi = phi;
BE_M_method.alpha = alpha;
BE_M_method.Re = Re;
BE_M_method.aa = wa;
BE_M_method.at = wt;
%save BE_M_method BE_M_method
