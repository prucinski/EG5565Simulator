%repeat (or start) calculation, status unknown
setAppToStatus(-1, app);

%% Tank dimensions, as per mechanical team, in meters
%default tank
if(app.NewTankCheckBox == false)  
    height = 2.4;
    width = 2.2;
    lengthOfElipsoidOnEnd = 0.4;
    lengthOfCylindrical = 5.858 - 2*lengthOfElipsoidOnEnd;
%custom tank
else
    height = app.HeightEditField.Value;
    width = app.WidthEditField.Value;
    lengthOfElipsoidOnEnd = app.EllipsoidEndRadiusEditField.Value;
    lengthOfCylindrical = app.TotallLengthEditField.Value - 2*lengthOfElipsoidOnEnd;
end
%% Sensor design dimensions, as per electrical team - kind of TODO
spacing = 7.75 * 1e-3; % mm * m
n = 1e-4;
braggWavelength = getBraggWavelength(spacing, n);   %1550 with arbitrary numbers
manufacturingTemp = 20;
%add wavelength shift to the graph
wavelengthShift = braggWavelength - app.peakWavelength;     
wavelengthShiftRef = braggWavelength - app.peakRefWavelength;
%display in the graph
textPosition = [0.95, 0.9]; % Adjust text position as needed
text(app.UIAxes, textPosition(1), textPosition(2), sprintf('Shift (temp/temp+strain): %.2f/%.2f nm', wavelengthShift*1e9, wavelengthShiftRef*1e9), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'right');
%% determine total volume of the tank first
%ellipse surface area: minor axis * major axis * pi
volumeOfCylindricalSection = height/2 * width/2 * pi * lengthOfCylindrical;
%On ends, there are two halfs of an elipsoid that together actually
%make a full elipsoid
volumeOfElipsoidSection = 4/3* pi *height/2 * width/2 * lengthOfElipsoidOnEnd;
totalVolume = volumeOfCylindricalSection + volumeOfElipsoidSection;
totalVolume = totalVolume * 1000; %conversion from m^3 to liters

%% Calculate the height of liquid based off of a height/volume input

%TODO: Height is meant to be acquired from spectral data

%Option 1 - liquid height
if(strcmp(app.ChooseInputDropDown.Value, app.ChooseInputDropDown.Items(1)))
    liquidHeight = app.ValueEditField.Value/1000;
    newVolume =  getEllipsoidPartialVolume(liquidHeight, height, width, lengthOfCylindrical, lengthOfElipsoidOnEnd);
    
%Option 2 - liquid volume
elseif(strcmp(app.ChooseInputDropDown.Value, app.ChooseInputDropDown.Items(2)))
    newVolume = app.ValueEditField.Value;
end
disp(['total volume of Tank:', num2str(totalVolume), ...
    ', current volume of Tank', num2str(newVolume)]); %apparently it's around 23186.7 liters
%paste the new volume into the application
app.VolumeGauge.Value = newVolume/totalVolume * 100;
app.VolumelitersEditField.Value = newVolume;
if(app.VolumeGauge.Value > 100)
   setAppToStatus(0, app); 
   errordlg('Volume levels above 100%','Volume Error');
end

%% Calculate the mass of the liquid
density = app.SpecificGravityEditField.Value;
mass = density * newVolume;
app.TotalMasskgEditField.Value = mass;



%% Finally, retrieve information back from the two spectra









