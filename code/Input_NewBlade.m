%WORLD CONSTANTS
WorldConstants

%CONDITIONS
Nr = inp.Nr;                   % number of blade elements

%GEOMETRY
if inp.normalRun disp(['Loading blade file ' inp.file_blade]); end
figure(100); text(0,0,{['Blade file:'],[strrep(inp.file_blade,'_','\_')]},'fontSize',28);
axis([-0.1 0.5 -0.2 0.5])
axis off
print('-dpdf',inp.name,'-append')
[blade] = readbladeS(inp); % reading the blade geometery from inp.file_blade
%disp('using 2*pi*R/(r - r0)*tan(beta0) distribution')
%Dm = 0.3; R = 0.15
%blade.r = linspace(R/2,R,Nr);
%blade.c = R-r;
r = blade.r;
chord = blade.c;
beta = blade.beta;
Dm = blade.Dm;
dr = blade.dr;
Nb = blade.Nb;

loc = isnan(chord); chord(loc) = 0; beta(loc) = 0;
R = Dm/2; r0 = min(r);

% positioning angle fault (or just difference), for entire blade
beta0 = inp.beta0;
blade.beta = (beta+beta0)*pi/180; %[rad]
