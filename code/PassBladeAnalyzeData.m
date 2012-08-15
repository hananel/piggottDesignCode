function inp = PassBladeAnalyzeData(inp)

% assuming tempCt.csv and tempCp.csv are outputs - 
% 2D matrices in VfVec and TSRvec of Ct and Cp respectively

CtData = csvread('tempCt.csv');
CpData = csvread('tempCp.csv');

inp.Ct_filename = input('insert CtData file name: ');
inp.Cp_filename = input('insert CpData file name: ');

% pointing to the CpCt directory
inp.CpCt_fileDirectory = [inp.file_input_directory 'CpCt/'];
csvwrite([inp.CpCt_fileDirectory '/' inp.Ct_filename], CtData);
csvwrite([inp.CpCt_fileDirectory '/' inp.Cp_filename], CpData);

% screen out
disp(['Cp and Ct data moved to ' inp.CpCt_fileDirectory ' under the file names ' inp.Ct_filename ' and ' inp.Cp_filename])