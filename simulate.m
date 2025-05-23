function [] = simulate(app)
    %repeat (or start) calculation
    %status unknown
    setAppToStatus(-1, app);
    %% Tank dimensions, as per mechanical team, in meters
    %default tank
    if(app.NewTankCheckBox.Value == false)  
        height = 2.6;
        width = 2.4;
        lengthOfElipsoidOnEnd = 1.05;
        lengthOfCylindrical =6.058 - 2*lengthOfElipsoidOnEnd;
        maxStress = 250;
    %custom tank
    else
        height = app.HeightEditField.Value;
        width = app.WidthEditField.Value;
        lengthOfElipsoidOnEnd = app.EllipsoidEndRadiusEditField.Value;
        lengthOfCylindrical = app.TotallLengthEditField.Value - 2*lengthOfElipsoidOnEnd;
        maxStress = app.MaxStressEditField.Value;
    end
    %% Sensor design dimensions, as per the user input (defaults per our design)
    [E_f, E_h, r_f, r_p, b_rp, h, a_f, a_h, G_p, G_a, n_eff, g_period,to_c,ref_T, v_f, L_f] = getTermsFromTable(app);
    %if using custom parameters:
    if(app.CustomHSMparametersCheckBox.Value == true)
        spacing = g_period*1e-6;         %convert from um to m; 
        n = n_eff;
        to_c = to_c*1e-6;  %convert from ue/K to e/K
        a_f = a_f *1e-6;    %ue/K to e/K

    %otherwise, we're using the same parameters as in the simulation
    else
        spacing = app.GratingPeriodEditField.Value;
        n = app.RefIndexEditField.Value;
        to_c =app.ThermoOpticEditField.Value;  
        a_f = app.ThermalCoeffEditField.Value; 
        L_f = app.GratingLengthEditField.Value;
    end
    braggWavelength = getBraggWavelength(spacing, n);   %1550 with default design
    %add wavelength shift to the graph
    wavelengthShift =  app.peakWavelength - braggWavelength;        %app.peakWavelength is inferred when loading in the file     
    wavelengthShiftRef = app.peakRefWavelength - braggWavelength;
    %display in the graph
    textPosition = [0.95, 0.9]; % Adjust text position as needed
    text(app.UIAxes, textPosition(1), textPosition(2), sprintf('Shift (temp/temp+strain): %.2f/%.2f nm', wavelengthShiftRef*1e9, wavelengthShift*1e9), ...
            'Units', 'Normalized', 'HorizontalAlignment', 'right');
    %% determine total volume of the tank first
    %ellipse surface area: minor axis * major axis * pi
    volumeOfCylindricalSection = height/2 * width/2 * pi * lengthOfCylindrical;
    %On ends, there are two halfs of an elipsoid that together actually
    %make a full elipsoid
    volumeOfElipsoidSection = 4/3* pi *height/2 * width/2 * lengthOfElipsoidOnEnd;
    totalVolume = volumeOfCylindricalSection + volumeOfElipsoidSection;
    totalVolume = totalVolume * 1000; %conversion from m^3 to liters
      
    %% Then, retrieve information back from the two spectra   
    p11 = 0.121;             %Pockel's constant (for calculation of k_e, k_t)
    p12 = 0.27;              %Pockel's constant
  
    k_e = 1 - 0.5*n^2*((1-v_f)*p12 - v_f*p11); 
    k_t = (1 - 0.5*n^2*(p11 + 2*p12))*a_f + to_c/n;

    dT = wavelengthShiftRef/(k_t*braggWavelength);
    app.TemperatureEditField.Value = ref_T + dT;    %reference spectrum is easy

   %Now, retrieve information from the strained spectrum
   %if simple model is chosen, recovering is easy
    if(app.AdvancedFBGmodelCheckBox.Value == false)
        strain = (wavelengthShift - k_t*braggWavelength*dT)/(k_e*braggWavelength);
        app.FBGStrainEditField.Value = strain;
        app.TankstrainEditField.Value = strain;
    %advanced model is trickier - for reference see transferMatrixMethod.m
    else
        [term0, lambdaTerm] = getLambdaTerm(E_f,E_h,G_p,G_a,r_f,r_p,b_rp,h);
      
       
        %recover e_t. a_f turned back to ue/K (requirement of the function)
        %strain passed doesn't matter so 0.01
        [e_t, ~] = getThermalAndMechanical(dT, 0.01, 2*L_f/3, 0, term0, lambdaTerm, E_f, a_h, a_f*1e6);
        %disp(e_t);
        %recover e_f (e_m_recovered) now
        e_m_recovered = (wavelengthShift - (k_t*dT + k_e*e_t)*braggWavelength)/(k_e*braggWavelength);
        app.FBGStrainEditField.Value = e_m_recovered;
      
        % cosh(lambdaTerm*0) = 1. Strongest component of the spectrum is at
        % it's strongest-bonded part (in the middle). E_f converted to Pa
        e_actual = e_m_recovered*E_f*1e9*term0/(1 - 1/cosh(lambdaTerm*2*L_f/3)); 
        app.TankstrainEditField.Value = e_actual;
        strain = e_actual;
    end



    %% Calculate the height of liquid based off of strain recovered
    % as per data provided by the Mechanical team   
    a = 5474.2;
    b = -0.0118;
    ratio = 1.2248/app.SpecificGravityEditField.Value;

    liquidHeight = ratio*(a*strain + b);
    disp(liquidHeight);
    %disp(height);
    if(liquidHeight > height)
        errordlg(["Warning: current strain readings indicate the liquid levels are above 100%, showing a height of", num2str(liquidHeight)], "Warning");
        liquidHeight = height;
    elseif( liquidHeight < 0)
        errordlg(["Warning: current strain readings indicate the liquid levels are below 0%, showing a height of", num2str(liquidHeight)], "Warning");
        liquidHeight = 0;
    end
    %disp(liquidHeight);
    newVolume =  getEllipsoidPartialVolume(liquidHeight, height, width, lengthOfCylindrical, lengthOfElipsoidOnEnd);
    
    %disp(['total volume of Tank:', num2str(totalVolume), ...
    %    ', current volume of Tank', num2str(newVolume)]); 
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

    %% Calculate maximum stress experienced by tank and FoS
    %see if stacked or not to see which formula to use
    %if(app.TankstackingCheckBox.Value == true)
    %    m = 799777027673.439;
    %    c = 13499481.179;
    %else %not stacked
        m = 608170824445.238;
        c = 1283290.86918642;
    %end
    stress = m*strain + c;
    stress = stress/1e6;  %convert Pa to MPa
    app.StressMPaGauge.Value = stress;
    app.StressMPaEditField.Value = stress;
    app.FactorofSafetyEditField.Value = maxStress/stress;

    %% Finally, set the app to status depending on parameters.
    if(app.FactorofSafetyEditField.Value < 1)
        setAppToStatus(0, app);
    elseif(app.FactorofSafetyEditField.Value < 1.2 || app.VolumeGauge.Value <75)
        setAppToStatus(1, app);
    else
        setAppToStatus(2, app);
    end










