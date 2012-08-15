function [P,RPM] = PMG(Vb,Kg,Rg,Rw,a,b)
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

RPMcut = Vb/Kg; % Vb = 24; Kg = 24/100  
Powerslope = Kg*Vb/Rtot;% Rr = 0.288

if nargin == 6
    if a>RPMcut
        RPM = a:100:b;
    else
        RPM = [a RPMcut:100:b];
    end
    P = Powerslope.*(RPM-RPMcut);
    for i=1:length(P)
        if P(i)<0 P(i)=0;
        end
    end
else
    P = Powerslope*(a - RPMcut); %[watt]
    if RPM<RPMcut
        P=0;
    end
end
