function  BE_M_method =  ActuatorDiscH(inp,blade)
% BEM actuator disk model.
% complementry to NewBlade.m
% uses NACA4412.m to calculate Cl and Cd as a function of alpha, according to table
% no reynolds related changes where accounted for at this time
% 30.8.10

Vf = inp.Vf;
TSR = inp.TSR;
draw = inp.doplots;
iter = inp.iter;
maxiter = inp.maxiter;
relax = inp.relax;

Nb = blade.Nb;
R = blade.Dm/2;
r = blade.r;
dr = blade.dr;
beta = blade.beta;
chord = blade.c;

%WORLD CONSTANTS
WorldConstants

% turbine properties
Nr = length(r);              % number of blade elements
vtip = Vf*TSR;
W  = vtip/(R);               % [radv/sec] rotational speed
sigma = chord.*Nb./(2*pi.*r);
ac = 0.2;
% initializing averge variables, for hysterysis in the convergence due to
% jump in 2D section data
avgcount = 1;
avgdT = 0;
avgdP = 0;
avgphi = 0;
avgalpha = 0;
avgRe = 0;
avgaa = 0;
avgat = 0;

if isfield(inp,'aa')
    aa = real(inp.aa);
    at = real(inp.at);
else
    aa = 0*ones(Nr,1);
    at = 0*Vf*r.*ones(Nr,1);
end

% convergence criteria
epsilon = inp.epsilon;
iter_criteria = 1;
convergence_criteria_aa = 100;
convergence_criteria_at = 100;
while and(or(max(convergence_criteria_aa)>epsilon,max(convergence_criteria_at)>epsilon),iter_criteria)
    
    fprintf(1,'.'); %fflush(stdout);
    
    phi = atan((1-aa)*Vf./((1+at).*W.*r));
    sinPhi = sin(phi);
    cosPhi = cos(phi);
    alpha = phi - beta;
    alphaDeg = alpha / pi*180;
    Vbe = Vf.*(1-aa)./sinPhi;
    Re = rho/ni.*Vbe.*chord;

    Re = killnans(Re,1); alphaDeg = killnans(alphaDeg,1);
    
    % calculating blade element properties
    [Cl Cd] = LoadSectionS(real(alphaDeg), real(Re), blade); 

    Cn = Cl.*cosPhi+Cd.*sinPhi;
    Ct = Cl.*sinPhi-Cd.*cosPhi;
    
    f = Nb/2.*(R-r)./(r.*sinPhi);
    F = 2/pi*acos(exp(-f));
    
    % correcting for aa>ac
    aa_new = 1./(4.*F.*sinPhi.^2./(sigma.*Cn) + 1);
    for j=1:Nr
        if aa_new (j)>ac
            K = 4*F(j).*sinPhi(j).^2./(sigma(j).*Cn(j));
            aa_new(j) = 0.5*(2+K.*(1-2*ac)-sqrt((K*(1-2*ac)+2).^2+4*(K*ac^2-1)));
        end
    end
    at_new = 1./(4.*F.*sinPhi.*cosPhi./(sigma.*Ct) - 1);
    
    %relaxation and sustitution
    aa = aa_new + relax*(aa - aa_new);
    at = at_new + relax*(at - at_new);
    
    %convergence
    convergence_criteria_aa = abs(max((aa_new-aa)./aa));
    convergence_criteria_at = abs(max((at_new-at)./at));
    
    if iter>maxiter
        iter_criteria=0;
    end
    iter = iter+1;
    
    % forces
    %Pn = Cn.*0.5*rho.*Vbe.^2.*chord;
    %Pt = Ct.*0.5*rho.*Vbe.^2.*chord;
    
    dT = 4*pi*r.*rho.*Vf^2.*aa.*(1-aa).*F*dr;
    dP = 4*pi*r.^3.*rho.*Vf.*W^2.*at.*(1-aa).*F*dr;
    if draw
        ff = figure(1); set(ff,'Position',[50 100 600 500 ])
        subplot(3,2,3);
        title('red: a_t , black: a\_a ');
        plot(r,at,'r',r,aa,'k');hold off
        axis tight;
        
        subplot(3,2,[1 2]);
        plot(r,dP,'-r',r,dT,'k');legend('dP','dT','Location','Best');
        %axis([0, r(end), min(dCp) ,max(dCp)]);
        title(['Moments and forces along blade. || TSR = ' num2str(TSR) ' V = ' num2str(Vf) ' m/s iteration #' num2str(iter)]);
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

% force calculation
%A = (Pt(2:end)-Pt(1:end-1)) ./ (r(2:end)-r(1:end-1));
%B = (Pt(1:end-1).*r(2:end)-Pt(2:end).*r(1:end-1)) ./ (r(2:end)-r(1:end-1));
%dP = Nb.*W*(1/3*A.*(r(2:end).^3-r(1:end-1).^3) + 1/2*B.*(r(2:end).^2-r(1:end-1).^2));

%A = (Pn(2:end)-Pn(1:end-1)) ./ (r(2:end)-r(1:end-1));
%B = (Pn(1:end-1).*r(2:end)-Pn(2:end).*r(1:end-1)) ./ (r(2:end)-r(1:end-1));
%dT = Nb*(1/3*A.*(r(2:end).^3-r(1:end-1).^3) + 1/2*B.*(r(2:end).^2-r(1:end-1).^2));

% dT = 0.5*rho.*Vbe.^2*Nb.*chord.*(Cl.*cosPhi+Cd.*sinPhi)*dr;              %eq. (7)
% dP = 0.5*rho.*Vbe.^2*Nb.*chord.*(Cl.*sinPhi-Cd.*cosPhi)*W.*r*dr;         %eq. (8)
%err_at = at - 1./(4*pi*rho.*r.^3*Vf*W^2.*F.*(1-F.*aa)).*(dP./dr)
%err_aa = aa - (1-sqrt(1-dT./dr./(pi*rho.*r*Vf^2)))./(2.*F)

if sum(isnan(dP))>0
    dP = killNans(dP);
    dT = killNans(dT);
end
BE_M_method.dP = dP;
BE_M_method.dT = dT;
BE_M_method.phi = phi;
BE_M_method.alpha = alpha;
BE_M_method.Re = Re;
BE_M_method.aa = aa;
BE_M_method.at = at;
