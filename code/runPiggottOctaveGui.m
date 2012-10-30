function returnValue = runPiggottOctaveGui(flag)
% using the zenity package, the simple gui commands are used as a documentation of this code
% install the zenity pkg from http://octave.sourceforge.net/zenity/ by downloading, opening octave and running pkg install zenity* (where * stands for the version name). This was written with 0.5.7
% run without flag
%
% written by Hanan Einav-Levy 17-08-2012

clear all
close all
warning off

pkg load zenity 

if nargin<1 
 	zenity_text_info("Piggott turbine design code", \
				 	["Welcome to the simple Octave Gui for piggott-turbine-design, a git hub project (tobe 17-08-2012)\n" \
				 	"Run this gui from the 'code' directory where all m files and inp_*.m files are located\n" \
				 	"The gui will go through the design stages of a Hugh Piggott HAWT. For detailes on the design see http://scoraigwind.co.uk/\n" \
				 	"\n" \
				 	"silly remark - double click on list items, other wise zenity crashes (do not use ok button\n" \
				 	"\n" \
				 	"1. First stage - design of blades. choose your inp_*.m file from the next list. Notice:" \
				 	" if you want to create a new input file, copy an inp_*.m file and keep the same inp_*.m name"]);
end

inpFile = zenity_list("select input file", {"Name"},{dir('inp_*.m').name});

% load input
run(inpFile);
figure(100); 
text(0,0.4,datestr(now),'color','r','fontSize',24)
text(0,0.2,{['Input file:'],strrep(inpFile,'_','\_')},'fontSize',28)

runBladeElement = zenity_list("run blade element?",{"option"},{"yes","no"});

% run blade element?
switch runBladeElement
 	case {"yes" }
    	returnValue = NewBlade(inp);         
	case {"no"}
		isFile = dir('../Output/Cp_CometME4.2_Nb3_NACA44XX.csv');
		if isempty(isFile)
			runBladeElement = zenity_list("you don't have a Cp or Ct file",{"your new option"},{"run-blade-element","choose-another-inp-file"},"all");
			switch runBladeElement
 				case {"run-blade-element" }
    				[minimum] = NewBlade(inp);         
				case {"choose-another-inp-file"}
					returnValue = feval("runPiggottOctaveGui",1)
			endswitch
		end
endswitch

% run power curve?
runPowerCurve = zenity_list("2. match generator to blade Cp curve?",{"option"},{"yes","no,I-wish-to-redesign-the-blades."});
switch runPowerCurve
 	case {"yes" }
    	returnValue = PowerCurve(inp);        
	otherwise
		returnValue = feval("runPiggottOctaveGui",1)
endswitch
