%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blade analyze input files
inp.file_input_directory = '/Users/hananlevy/src/piggott-turbine-design/octave/Input/';
inp.file_Output_directory = '/Users/hananlevy/src/piggott-turbine-design/octave/Output/';
inp.file_blade_directory = [inp.file_input_directory 'Geometry/'];
inp.file_blade = 'CometME4.2_Nb3_NACA44XX.txt';%'WindyLightGO407A_Nb4_shavshevet_TSR3_clean_5_0.txt'; %'Thomas-2,1m_Nb3.csv';
% 'WindyLightGO407A_Nb4_shavshevet_TSR3_clean_10_0_r0_0_06_R.txt'; 
%'CometME25_Nb3.txt'; %'WindyLightGO407A_Nb4_shavshevet_TSR3.txt'; 
%'CometME25_Nb3.txt'; %
% blade analyze physical parameters
inp.beta0 = 0;
inp.Cd0 = 0; % 0.05 windylight last
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
inp.Vfend = 10;
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

%diode option - Adital Ela
inp.diode = 0;
inp.diode_p = 4;
inp.diode_s = 2;
inp.I0 = 0.2147 * inp.diode_p;  
inp.t0 = 0.8451 / inp.diode_s;  
inp.v0 = 2.4470 * inp.diode_s;

% Test data
inp.compare2TestData = 1;
inp.TestData_fileDirectory = [inp.file_input_directory 'TestData/'];
inp.TestData_file = 'CometME4.2_eff.csv'; %'CometME3_data.txt';% 'WA2_Sep2010.csv'; %

inp.plotoff = 0;
inp.chord75 = 0.08;             %[m]
 
% generator data
inp.coil.series = 5;            % number of coils at series
inp.coil.N = 15;                 % number of coils
inp.coil.FF = 0.62;             % actual filling of coil rectangle with wire. Hugh Piggot's number
inp.coil.turns = 54;
inp.coil.thickness = 12;
inp.magnet.N = 20; %12          % number of magnets
inp.magnet.width = 30;          % [mm]
inp.magnet.length = 46;         % [mm]
inp.magnet.B = 0.72;            % [Tesla] - Hugh uses 0.62 for two disk configuration. 
inp.Rstator = 0.202;            % [m] Radius to middle of magnets/coils Noam 3meter - 0.152

% mesured data for the PMG
inp.measuredPMG = 1;
inp.Kg = 0.436; % 1/100; %                       % 285791 p. 80     % 1/100;       % [Volt / RPM]
inp.Rg = 1.58; %7.46; %                        % 7.46 [Ohm] 
% motor conventions
% Ke = inp.Kg back EMF constant (normally in mV/RPM)
% Km = torque constant in mNm/A
% Km = Ke*30/pi
% 
% wiring and battery data
inp.Vb = 48; %24                                        % [Volt]
% load controlVoltage 
inp.useVoltageControl = 0; % 0 = don't use, 1 = usee optimal power control, 2 use optimal TSR control
inp.Rw = 0.005;%0.01;             %  [Ohm]
inp.r = 0.05;%0.05;             % [Ohm] battery internal resistence

%looping through relevant wind speeds
inp.VvecStep = 1;                                        %[m/s]

%Color for plots
inp.color = 'b';

% plot options
inp.plotGen = 1;
if inp.measuredPMG
    inp.plotGen = 0;
end

disp('input taken from inp_default')

% windylight solution maxon motors p.80 90 watt motor number 285791
% Kg = 1/100 volt/RPM, Rg = 7.46 ohm

