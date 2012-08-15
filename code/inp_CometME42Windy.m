%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blade analyze input files
inp.normalRun = 1;

% Blade analyze input files
inp.file_input_directory = '/home/hanan/Documents/src/piggott-turbine-design/octave/Input/';
inp.file_Output_directory = '/home/hanan/src/piggott-turbine-design/octave/Output/';
inp.file_blade_directory = [inp.file_input_directory 'Geometry/'];
inp.file_blade = 'CometME4.2_Nb3_NACA44XX.txt';

% blade analyze physical parameters
inp.beta0 = 0;
inp.Cd0 = 0; 

% Blade analyze model
inp.useOhad = 1;

% Blade analyze numerical parameters
inp.Nr = 50;
inp.param = 0;
inp.iter = 1;
inp.maxiter = 40;
inp.relax = 0.5;
inp.epsilon = 0.01;

% Blade analyze range
inp.TSRstart = 1; 
inp.TSRend = 14; 
inp.TSRstep = 1;
inp.Vfstart = 2;
inp.Vfstep = 1;
inp.Vfend = 9;

% Blade analyze output parameters and file names
inp.filename = 'temp';
inp.doplots = 0;
inp.doplotsCp = 1;
inp.CpCt_fileDirectory = 'Output/';
inp.CpCt_filename = ['CpCt_' strrep(inp.file_blade,'.txt','')]; % graphic file
inp.Ct_filename =  ['Ct_' strrep(inp.file_blade,'txt','csv')];
inp.Cp_filename = ['Cp_' strrep(inp.file_blade,'txt','csv')];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power Curve analyze Input

%diode option - Adital Ela windyLight only
inp.diode = 0;
inp.diode_p = 4;
inp.diode_s = 2;
inp.I0 = 0.2147 * inp.diode_p;  
inp.t0 = 0.8451 / inp.diode_s;  
inp.v0 = 2.4470 * inp.diode_s;

% Test data
inp.compare2TestData = 1;
inp.TestData_fileDirectory = [inp.file_input_directory 'TestData/'];
inp.TestData_file = 'power_wind_shaeb_windyBoy.txt';
inp.R = 2.1; % [m]

% plotting
inp.plotoff = 0;

% misc
inp.chord75 = 0.08;             %[m]
 
% generator data
inp.coil.series = 5;            % number of coils at series
inp.coil.N = 15;                 % number of coils
inp.coil.FF = 0.75;             % actual filling of coil rectangle with wire. Hugh Piggot's number
inp.coil.turns = 240;
inp.coil.thickness = 12;
inp.magnet.N = 20; %12          % number of magnets
inp.magnet.width = 30;          % [mm]
inp.magnet.length = 46;         % [mm]
inp.magnet.B = 0.72;            % [Tesla] - Hugh uses 0.62 for two disk configuration. 
inp.Rstator = 0.202;            % [m] Radius to middle of magnets/coils Noam 3meter - 0.152
inp.coil.factor = 1.3; %impedence factor - empirical correction according to Hugh Piggott's book

% mesured data for the PMG
inp.measuredPMG = 1;
inp.Kg = 1.7633; % [Volt / RPM]
inp.Rg = 20.8 ;  % [Ohm] 
inp.RPMcor = 1;  % a empirical correction for reduced Kg with RPM - multiplying RPM by linspace(1,RPMcor,N) where N=length(RPM)
% motor conventions
% Ke = inp.Kg back EMF constant (normally in mV/RPM)
% Km = torque constant in mNm/A
% Km = Ke*30/pi
% 
% wiring and battery data
inp.Vb = 575.6669; %24                                        % [Volt]
inp.loadVlaw = 0; % law of voltage change with current - using proxy for velocity instead, 
                  % based on simple piece wise linear interpolation from CometME4.2 data
% load controlVoltage 
inp.WindyBoy = 1;
inp.useVoltageControl = 1; % 0 = don't use, 1 = usee optimal power control, 2 use optimal TSR control
inp.VoltageControlDelta = 0;
inp.Rw = 0.05;            %  [Ohm]
inp.r = 0;             % [Ohm] battery internal resistence

%looping through relevant wind speeds
inp.VvecStep = 1;                                        %[m/s]

%Color for plots
inp.color = 'g';

% plot options
inp.plotGen = 0;
if inp.measuredPMG
    inp.plotGen = 0;
end

disp('input taken from inp_CometME42Windy')

% windyBoy control law
% [m/s] [Watt]      [Volt]
inp.controlSpeed = [2:11 20];
inp.controlVolt =  [ 104.7600  154.8444  202.7304  249.1157  293.3207  336.0024  376.6624  431.4216 490.1566 553.0192 553.0192];
%inp.controlVolt =  [  106.2445  155.8352  202.7124  247.2475  289.2770  346.2175 407.3966  470.8493  533.7673  575.6669 575.6669];
%ControlLaw1 = [  2     21.68       104.76
%                3     72.0817     154.8444
%                4     167.4486    202.7304
%                5     321.539     249.1157 
%                6     544.7512    293.3207
%                7     849.6605 	336.0024
%                8     1243.966 	376.6624
%                9     1737.349 	431.4216
%                10    2340.5397   490.1566
%                10.9  2986.7802   553.0192
%                20      nan       553.0192];
ControlLaw2 = [ 2 20.3188 106.2445
                3 66.6575 155.8352
                4 153.5471 202.7124
                5 292.2203 247.2475
                6 491.943 289.277
                7 761.104 346.2175
                8 1109.3761 407.3966
                9 1545.5199 470.8493
                10 2077.845 533.7673
                11 2711.3249 575.6669];
% sensitivity test for control law
DeltaV = [-0.2 -0.1 0 0.1 0.2];                            % difference in control law - this is the precent change in voltage set at each wind speed (or power) 
DeltaE = [0.9555    0.9891    1.0000    0.9864    0.9624]; % difference in energy for a k=2 A=5 (Vavg~6 m/s) site
% also, different between controlLaw1 and 2 is negligible

% remarks about measurements 14/4/12
% 4.2 m	18 m tower	wind sensor at 15 m
% V(DC)/RPM=2 @open circuit		
% 3 phase		
% 5 coils per phase		
% 240 turns per coil		
% Data is an averaged of 15 m		

save InpFile inp
