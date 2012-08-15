%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blade analyze input files
inp.normalRun = 1;

% Input and output file locations
inp.file_input_directory = '/Users/hananlevy/src/piggott-turbine-design/octave/Input/';
inp.file_Output_directory = '/Users/hananlevy/src/piggott-turbine-design/octave/Output/';
inp.file_blade_directory = [inp.file_input_directory 'Geometry/'];
inp.file_blade = 'CometME4.2_Nb3_NACA44XX.txt';
inp.R = 2.1; % [m]

% blade analyze physical parameters
inp.beta0 = 0;  % blade positioning angle (added to beta distribution)
inp.Cd0 = 0;    % constant drag of blade section - added to Cd(alpha). keep 0 normally.

% Blade analyze model
inp.useOhad = 1; % 1 for Ohad Gur 2009 BEM model, 0 for Hansen 2000 BEM model

% Blade analyze numerical parameters
inp.Nr = 10;        % number of blade elements
inp.param = 0;      % for the optimization rutine
inp.iter = 1;       % initial iteration number
inp.maxiter = 20;   % maximum iterations
inp.relax = 0.5;    % relaxation parameter
inp.epsilon = 0.01; % convergence criterion

% Blade analyze range
inp.TSRstart = 1; 
inp.TSRend = 14; 
inp.TSRstep = 1;
inp.Vfstart = 2;
inp.Vfstep = 1;
inp.Vfend = 12;

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

% diode option - for a diode load (no battery)
inp.diode = 0;
inp.diode_p = 4;
inp.diode_s = 2;
inp.I0 = 0.2147 * inp.diode_p;  
inp.t0 = 0.8451 / inp.diode_s;  
inp.v0 = 2.4470 * inp.diode_s;

% Test data
inp.compare2TestData = 1;
inp.TestData_fileDirectory = [inp.file_input_directory 'TestData/'];
inp.TestData_file = 'CometME4.2_P.txt'; %'CometME3_data.txt';% 'WA2_Sep2010.csv'; %

inp.plotoff = 0;
inp.chord75 = 0.08;             %[m]
 
% generator data
inp.coil.series = 5;            % number of coils at series
inp.coil.N = 15;                 % number of coils
inp.coil.FF = 0.86;             % actual filling of coil rectangle with wire. Hugh Piggot's number
inp.coil.turns = 70;
inp.coil.thickness = 12;
inp.magnet.N = 20; %12          % number of magnets
inp.magnet.width = 30;          % [mm]
inp.magnet.length = 46;         % [mm]
inp.magnet.B = 0.72;            % [Tesla] - Hugh uses 0.62 for two disk configuration. 
inp.Rstator = 0.202;            % [m] Radius to middle of magnets/coils Noam 3meter - 0.152
inp.coil.factor = 1.3; %impedence factor - empirical correction according to Hugh Piggott's book

% mesured data for the PMG
inp.measuredPMG = 0;    % if 1 - use measured data, if 0 - calculate from above parameters
inp.Kg = 0.4;   % [Volt / RPM]
inp.Rg = 0.68;  % [Ohm] 
inp.RPMcor = 1; % a empirical correction for reduced Kg with RPM - multiplying RPM by linspace(1,RPMcor,N) where N=length(RPM)
% motor conventions
% Ke = inp.Kg back EMF constant (normally in mV/RPM)
% Km = torque constant in mNm/A
% Km = Ke*30/pi
% 
% wiring and battery data
inp.Vb = 48; %24         % [Volt]
inp.loadVlaw = 0; % law of voltage change with current - using proxy for velocity instead, 
                  % based on simple piece wise linear interpolation from CometME4.2 data
% load controlVoltage 
inp.WindyBoy = 0;
inp.useVoltageControl = 0; % 0 = don't use, 1 = usee optimal power control, 2 use optimal TSR control
inp.Rw = 0.01;%0.01;             %  [Ohm]
inp.r = 0.1;%0.05;             % [Ohm] battery internal resistence


%looping through relevant wind speeds
inp.VvecStep = 1;                                        %[m/s]

%Color for plots
inp.color = 'b';

% plot options
inp.plotGen = 1;
if inp.measuredPMG
    inp.plotGen = 0;
end
save InpFile inp
disp('input taken from inp_CometME42')

% windylight solution maxon motors p.80 90 watt motor number 285791
% Kg = 1/100 volt/RPM, Rg = 7.46 ohm

