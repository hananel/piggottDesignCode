function [Kg,Rg,coil] = PMGdesign(R,magnet,coil)
% according to Hugh Piggot's notes about stator design from the 2010 kindle
% edition of his axial windturbine plans + some geometrical considerations
% of axial stators
%
% assumptions - a 3 phase machine, and Kg is given for the average voltage
% R is the radius to the middle of the coils/magnets [ meters]
% magnet is a struct, containing
% magnet.N                               number of magnets
% magnet.width                    width of magnet [mm]
% magnet.length                  length of magnet [mm]
% magnet.B                               magnet field strength [Tesla] 
% the length of the magnet is assumed at the radial direction
% coil is a struct, that comes in to this function with the following information
% coil.N                                      total coil number
% coil.series                         amount of coils connected in series
% coil.turns                            number of turns in one coil
% coil.thickness                  thickness of coil
% coil.FF                                    coil form factor - efficiency
% of making it, volume of actual copper divided by volume of coil
%
% and leaves with the additional parameters:
% coil.width
% coil.legWidth
% coil.sqmm
% coil.lAvg
% coil.l
% coil.m
% coil.R
% coil.S
% coil.Pmax
% coil.cooling

coil.width            = 2*pi*(R-magnet.length/1000/2) / (coil.N)*1000;       %[mm]
coil.legWidth         = (coil.width - magnet.width)/2;                       %[mm]
coil.lAvg             = 2*(magnet.width+magnet.length) + pi * coil.legWidth; %[mm]
coil.l                = coil.lAvg * coil.turns;
coil.sqmm             = coil.legWidth* coil.thickness * coil.FF / coil.turns;
if coil.AWG
    [coil.AWG,coil.sqmm,coil.AWG_sqmm,coil.wireN]  = AWGfinder(coil.sqmm);
else
    if coil.sqmm>pi*(1.8/2^2) % Piggott's plans seem to use this rule
        coil.wireN = 2;
    else
        coil.wireN = 1;
    end
    coil.mm = 2*round(sqrt(coil.sqmm/coil.wireN/pi)*10)/10; % diameter [mm]
    coil.sqmm = coil.wireN*pi/4*coil.mm^2; % [mm^2]
end
coil.m                = coil.l * coil.sqmm * 0.009;                          %[gr]
coil.R                = coil.l / coil.sqmm / 56000;                          %[Ohm @ 20C]
Rg                    = 2 * coil.series * coil.R * coil.factor;      %[Ohm @ 20C] - corrected by 1.3 thumb rule, to imply impedence rather then resistence
Kg                    = magnet.width/1000 * magnet.length/1000 * magnet.N *  magnet.B * coil.series * coil.turns * 1.73 / 30; %[V/RPM] 
coil.Rg = Rg; coil.Kg = Kg;

% the 1.73 factor comes from the line voltage being  1.73 times the phase voltage since we have 3 phases
% In Hugh notes there is an additional multiplication by 1.57 that comes from: 
% peak voltage being 1.57 times the average voltage, and 
% Vbatt = 24;
% cutInRPM = (Vbatt+1.4)/Kg

% From battery to diode box three phases cable of 16 SQR mm for 12 meters, 
% from diode box to battery for each polarity two cables of 35 SQR mm (i.e 70 SQR mm for the (+) and 70 SQR MM for the (-). 
% Distance from rectifier box to battery is around 40 meters.
