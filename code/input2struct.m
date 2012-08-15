function [ inp ] = input2struct( input_type, input_args )
%INPUT2STRUCT takes input from web page and builds a matlab/octave struct
%   input_args are all input parameters included in the script
%   inp_deafualts.m [ and powerCurve_input.m ]
%   exceptions are input files, which could be saved and then the directory
%   and file will be transfered as other inputs, or could be changed to
%   just be the actual input - the vectors from the files

if input_type==0 % blade geometrey analysis
    inp.Nr =                     input_args(1);
    inp.useNACA =                input_args(2);
    inp.beta0 =                  input_args(3);
    inp.file_blade_directory =   input_args(4);
    inp.file_blade =             input_args(5);
    inp.doplots =                input_args(6);
    inp.doplotsCp =              input_args(7);
    inp.filename =               input_args(8);
    inp.TSRstart =               input_args(9);
    inp.TSRend =                 input_args(10);
    inp.TSRstep =                input_args(11);
    inp.Vf =                     input_args(12);
    inp.param =                  input_args(13);
    inp.useOhad =                input_args(14);
    inp.iter =                   input_args(15);
    inp.maxiter =                input_args(16);
    inp.relax =                  input_args(17);
    inp.epsilon =                input_args(18);
    inp.compare2TestData =       input_args(19);
    inp.TestData_file =          input_args(20);
    inp.TestData_fileDirectory = input_args(21);
    inp.rho =                    input_args(22); %[Kg/m^3]
    inp.R =                      input_args(23);
    inp.Cp_file =                input_args(24);
    inp.CpCt_fileDirectory =     input_args(25);
    inp.Ct_file =                input_args(26);
else            % turbine analysis
    inp.chord75 =                input_args(1);  % [m]      Chord of 3/4 radius for Reynolds calculation at each windspeed
    inp.coil.series =            input_args(1);  %          Number of coils at series
    inp.coil.N =                 input_args(2);  %          Number of coils
    inp.coil.FF =                input_args(3);  %          Actual filling of coil rectangle with wire. Hugh Piggot's number
    inp.coil.turns =             input_args(4);
    inp.coil.thickness =         input_args(5);
    inp.magnet.N =               input_args(6);  %          Number of magnets
    inp.magnet.width =           input_args(7);  % [mm]
    inp.magnet.length =          input_args(8);  % [mm]
    inp.magnet.B =               input_args(9);  % [Tesla]  Hugh uses 0.62 for two disk configuration.
    inp.Rstator =                input_args(10); % [m]      Radius to middle of magnets/coils
    inp.measuredPMG =            input_args(11); %          Mesured data for the PMG
    inp.Vb =                     input_args(12); % [Volt]   Battery volts
    inp.Rw =                     input_args(13); % [Ohm]    Wire resistence between generator and battery
    inp.r  =                     input_args(14); % [Ohm]    Battery internal resistence
    inp.Vvec =                   input_args(15); % [m/s]    Wind speed range
    inp.color =                  input_args(16); %          Color for plots
    inp.plotGen =                input_args(17); %          Plot options
    inp.useVoltageControl =      input_args(18); %          0 = don't use, 1 = usee optimal power control, 2 use optimal TSR control
    inp.plotoff =                input_args(19);
end

end

