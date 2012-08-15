function [NACA0012_shift] = shiftNACA0012_byFeature(Section, NACA0012)
% Section is a struct
% Section.alpha
% Section.Cl
% Section.Cd
% NACA0012 has the same cells, and also the outpout NACA0012_shift will
% have the same cells
% The measured post stall NACA0012 data is used for the alpha angles (AOA) that are higher then
% measured or xfoil simulated.
% This function shifts the NACA0012 data to more closely resemble the
% Section it is used with, so that no discontinuties arise (which interfere
% with convergence, and are probably less physical)
%
%
% The method used is shifting the Cl trend according to the stall point,
% just by the max points alighning, and then shifting up and down according
% to least square fitting around it just for the dCl.


% finding maximum point
[NACA_max,NACA_n] = max(NACA0012.cl);
[Cl_max,Cl_n] = max(Section.cl);
dCl = NACA_max - Cl_max;
da = NACA0012.alpha(NACA_n) - Section.alpha(Cl_n);
% shifting it
NACA0012_shift.alpha = NACA0012.alpha - da;
NACA0012_shift.cl = NACA0012.cl - dCl;
% TODO - minimum square fitting around this point

% TODO - smoothing the change between the end of Section Cl and
% NACA0012_shift

% TODO - shifting Cd as well by the same alpha
NACA0012_shift.cd = NACA0012.cd;

function fy = minSquarey(alpha, Cls, Cl0012, da, dCl)
fy = sqrt(sum(Cls.^2-(intepr1(alpha,Cl0012,alpha+da)+dCl).^2));
