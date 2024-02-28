%Calculate the change in strain based off of a model where strain and
%temperature are independent
%Method from here: https://iopscience.iop.org/article/10.1088/0957-0233/8/4/002

%Units: braggShift: um
%       dTemp: K
%       tempSensitivity: pm/K
%       strainSensitivity: pm/ue
%       where u = micro, e - strain
function dStrainSimple = getStrainFromSuperimposedFBG(braggShift, dTemp, tempSensitivity, strainSensitivity)
    braggShiftInPico = 1000000*braggShift;                                                      %consistent units
    dStrainSimple = (braggShiftInPico - dTemp*tempSensitivity)/strainSensitivity;               %in ue
end


%example from the paper: for 1.55 um, sS is 1.2 pm/ue, tS is 13 pm/C