function blade = readbladeS(inp,dummy)
% readbladeS is a general blade geometry reading function
% which includes reading the section parameters as well as the chort and
% twist. 
%%
% The function expects [inp.file_blade_directory '/' inp.file_blade]
%%
% to be a csv file (comma delaminated) containg 4 columns: 
% radius [m] chord [m] twist [deg] Section name "Name"
% for HughPiggot wooden blades, a preprocessing needs to be done to prepare
% the csv file, and according to the thickness to chord ratio (in precent -
% typically from 11-22, name the Sections as "NACA44XX" according to this
% thickness. for example - for thickness ratio of 15%, the section name
% should be "NACA4415"

% Read geometry file
if isstruct(inp)
    fileName = [inp.file_blade_directory inp.file_blade];
    Nr = inp.Nr;
else
    fileName = inp;
end

fid = fopen(fileName);
j = 1;
while 1
   lineRaw = fgetl(fid);
   if ~ischar(lineRaw), break, end
   lineOut = textscan(lineRaw,['%f %f %f %s'],'delimiter',',');
   blade.rOrig(j)      = cell2mat(lineOut(1));
   blade.c(j)   = cell2mat(lineOut(2));
   blade.beta(j)    = cell2mat(lineOut(3));
   blade.Section(j) = strrep(lineOut{4},'"','');
   j = j+1;
end

% process data and interpolate into blade elements

if nargin==2
    figure;
    [AX,H1,H2] = plotyy(blade.rOrig,blade.c,blade.rOrig,blade.beta);
    text(blade.rOrig,blade.c+max(blade.c)*0.1,blade.Section)
    ylabel(AX(1),'chord[m]'); ylabel(AX(2),'beta [deg]');
    xlabel('radius [m]');
    title(inp.file_blade)
    blade.r = blade.rOrig;  
    blade.dr = blade.r(2)-blade.r(1);
    blade.Dm = 2*blade.r(end);
    blade.Nb = str2num(inp.file_blade(strfind(inp.file_blade,'_Nb')+3));
else
    blade.Nb = str2num(inp.file_blade(strfind(inp.file_blade,'_Nb')+3));
    % derived properties
    R = blade.rOrig(end);                      % [m]  end of blade
    if isfield(inp,'r0')
        r0 = inp.r0;
    else
        r0 = blade.rOrig(1);                    % [m]  start of blade
    end
    dr = (R-r0)/Nr;
    blade.r = linspace(r0,R-dr,Nr)';            % radial coordinant of blade elemnts  -0.1*dr
    blade.c = interp1(blade.rOrig,blade.c,blade.r);
    blade.beta = interp1(blade.rOrig,blade.beta,blade.r);
    blade.dr = blade.r(2)-blade.r(1);
    blade.Dm = R*2;
end
fclose(fid);