%% Parameters
usingApp = -1;  %#ok<NASGU>   
lambdaTerm = -1; 
simpleMode = 1;
try %check if the script is being ran with app or without
    cla(app.FBGStrainedGraph);                         %clear the figure
    usingApp = 1;
catch
    usingApp = 0;
end
if(usingApp == 1)
    if(app.AdvancedFBGmodelCheckBox.Value == true)
        simpleMode = 0; %Is code using simple or advanced mode?
        [E_f, E_h, r_f, r_p, b_rp, h, a_f, a_h, G_p, G_a, ~, ~,~, ~, ~] = getTermsFromTable(app);
        [term0, lambdaTerm] = getLambdaTerm(E_f, E_h, G_p, G_a, r_f, r_p, b_rp, h);    %only needs calling once
    end
    core_refractive = app.RefIndexEditField.Value;     %just for reference of variables
    core_cladding = app.CladdingIndexEditField.Value;                               
    changeInRefractiveIndex = app.DRefIndexEditField.Value; 
    gratingPeriod = app.GratingPeriodEditField.Value;                 %m
    gratingLength = app.GratingLengthEditField.Value;                 %m 
    %photoelasticCoefficient = app.PeCoeffEditField.Value;
    poissonsRatio = app.PeCoeffEditField.Value;
    thermalExpansionCoefficient = app.ThermalCoeffEditField.Value;    %1/K
    thermoOpticCoefficient = app.ThermoOpticEditField.Value;          %1/K
    strainInSection = app.StrainEditField.Value;
    dTemp = app.DTempEditField.Value;

%for using without the app
else 
    core_refractive = 1.47;
    core_cladding = 1.457;                             
    changeInRefractiveIndex = 1e-4;
    gratingPeriod = 5.27212e-7;                   %m
    gratingLength = 10*1e-3;                      %m 
    %photoelasticCoefficient = 0.215; %the photoealastic coefficient can also be calculated from Eq 3 from
    %here: https://link.springer.com/article/10.1007/BF02323100
    poissonsRatio = 0.17; 
    thermalExpansionCoefficient = 0.5e-6;       %1/K
    thermoOpticCoefficient = 1.4e-5;           %1/K 
    strainInSection = 0.0005;
    dTemp = -5;

end

v_f = poissonsRatio;
%pe = photoelasticCoefficient;

n_1 = core_refractive;        %just for reference of variables
n_2 = core_cladding;  
dNeff = changeInRefractiveIndex;
gratingP = gratingPeriod;
L= gratingLength;  
braggL = getBraggWavelength(gratingP, n_1);         %m  
alpha = thermalExpansionCoefficient;

% other constants

p11 = 0.121;             %Pockel's constant (paper 4)
p12 = 0.27;              %Pockel's constant

k_e = 1 - 0.5*n_1^2*((1-v_f)*p12 - v_f*p11); 
k_t = (1 - 0.5*n_1^2*(p11 + 2*p12))*alpha + thermoOpticCoefficient/n_1;

disp(['K_e: ', sprintf('%.12f', k_e), ' k_t: ', sprintf('%.20f', k_t)]);

%maximum N
M = 2*n_1*gratingLength/braggL;

%paper1: https://ijcsi.org/papers/IJCSI-9-1-2-368-374.pdf
%paper2:https://www.sciencedirect.com/science/article/pii/S235271101630022X
%paper3: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=618377
%paper4: https://www.mdpi.com/1424-8220/20/15/4223 


%% Continue on

currentLarray = linspace(braggL-1e-9,braggL+1e-9,2000);
y_result = zeros(1,200);
j = 1;
%% reference spectrum
for currentL = currentLarray
    %% Equation for fiber propagation constant B and constants k and y (eq 12, paper1)
    B = 2*pi*n_1/currentL;
    dB = B - pi/gratingP;
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
    %% Prepare for the transfer matrix method
    N = 40;
    sectionL = L/N;
    
    T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
    T_12 = -1i*k/y*sinh(y*sectionL);
    T_21 = 1i*k/y*sinh(y*sectionL);
    T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
    
    T_i = [T_11, T_12;T_21, T_22];
    
    T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
    for i = 2:N
        T= T*T_i;
    end
    %% extract reflected & transmitted amplitude from boundary conditions
    % 4.8.3 & 4.8.4 from Fiber Bragg Gratings
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    shouldBeOne = real(transmissivity)^2+ imag(transmissivity)^2 + real(reflectivity)^2+ imag(reflectivity)^2;  %all OK, distance from origin is 1
    R = k^2*sinh(y*L)^2/(dB^2*sinh(y*L)^2 + k^2*cosh(y*L)^2);
    %y_result(end+1) = R;    this has weird arifacts for some reason
    y_result(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;


end
%% Plot the reference spectrum
if(usingApp)
    plot(app.FBGRefGraph, currentLarray, y_result);
    xlim(app.FBGRefGraph, [currentLarray(1) currentLarray(end)]);
    ylim(app.FBGRefGraph, "auto");
    textPosition = [0.95, 0.9]; % Adjust text position as needed
    text(app.FBGRefGraph, textPosition(1), textPosition(2), sprintf('Bragg wavelength: %.2f nm', braggL*1e9), 'Units', 'Normalized', 'HorizontalAlignment', 'right');
else
    plot(currentLarray, y_result);
end

%% Now, the strained spectrum 
y_result_strained = zeros(1,200);
j = 1;
strainTable = {};
for currentL = currentLarray 
    %% - change induced by strain & temperature
    %new grating period (eq 2, paper 2 + eq3, paper3)
    %Check if we're using a strain file
    if(usingApp)
        if(app.StressResponseOK.Enable == true)
            %load the file only once
            if(isempty(strainTable))   
                strainTable = readtable(app.stressPath + "\" + app.stressFile);
                N = height(strainTable);
                disp(N);
            end
        end
    else
      N = 40;
    end
    sectionL = (2/3)*L/N;
    tempstrains = zeros(N, 2);
    currentLengthFromBeginning = sectionL*N/2 - sectionL;           %for advanced mode
    for i = 2:N
        if(~isempty(strainTable))
           strainInSection = strainTable{i, 5};    %access the strain value 
        end    
        if(simpleMode == 1)
              newN = n_1; %note - there should be a change of the effective refractive index  
              newGratingP = gratingP*(1 + k_e*strainInSection + k_t*dTemp); 
               newBraggL = getBraggWavelength(newGratingP, newN);%note - this is only useful for uniform stress
        else %we're using the complex method to infer the spectrum
            if(i <= N/2) %we only need to do half of the spectrum, past the middle point it's symmetric
              [e_t, e_m] = getThermalAndMechanical(dTemp, strainInSection, sectionL*N, currentLengthFromBeginning,term0, lambdaTerm, E_f, a_h, a_f);
              tempstrains(i,:) = [e_t, e_m];

              newGratingP = gratingP*(1 + k_e*(e_t + e_m) + k_t*dTemp);
                   if(i == 2)
                    newN = n_1; %note - there should be a change of the effective refractive index  
                    newBraggL = getBraggWavelength(newGratingP, newN);%note - this is only useful for uniform stress
                  end
              currentLengthFromBeginning = currentLengthFromBeginning - sectionL;   %advance along
            else
              tempstrains(i,:) = tempstrains(end/2 + 1 - i + N/2,:);
              e_t = tempstrains(i,1); e_m = tempstrains(i,2);
              newGratingP = gratingP*(1 + k_e*(e_t + e_m) + k_t*dTemp);
            end
        end
        %
        
        
        B = 2*pi*core_refractive/currentL; 
        dB = B - pi/newGratingP;        
        k = pi*dNeff/currentL;
        y = sqrt(k^2 - dB^2);
        T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
        T_12 = -1i*k/y*sinh(y*sectionL);
        T_21 = 1i*k/y*sinh(y*sectionL);
        T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
        T_i = [T_11, T_12;T_21, T_22];
        if(i == 2)
            T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
            if(isempty(strainTable))
               T= T*T_i; %For uniform sections it doesn't matter
            end
        else  
            T= T*T_i;
        end
    end
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    y_result_strained(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;
end
disp(tempstrains);
newBraggL_strain = newBraggL;   %for plotting later
y_result_temp = zeros(1,200);
j = 1;
for currentL = currentLarray 
    %% Change induced by only the temperature on the same fibre (no mechanical strain)
    newGratingP = gratingP*( 1 + k_t*dTemp); 
    newN = n_1; %note - there should be a change of the effective refractive index
    newBraggL = getBraggWavelength(newGratingP, newN);%note - this is only useful for uniform stress
    B = 2*pi*core_refractive/currentL; 
    dB = B - pi/newGratingP;        
    k = pi*dNeff/currentL;
    y = sqrt(k^2 - dB^2);
    N = 20;
    sectionL = (1/3)*L/N;         %to distinguish between spectra, temperature spectrum will have half the length
    T_11 = cosh(y*sectionL) - 1i*dB/y*sinh(y*sectionL);
    T_12 = -1i*k/y*sinh(y*sectionL);
    T_21 = 1i*k/y*sinh(y*sectionL);
    T_22 = cosh(y*sectionL) + 1i*dB/y*sinh(y*sectionL);
    
    T_i = [T_11, T_12;T_21, T_22];
    T = T_i;  %begin with T_1, then multiply by T_2, then T_3... 
    for i = 1:N
        T= T*T_i;
    end
    transmissivity = 1/T(1,1);
    reflectivity = T(2,1)/T(1,1);
    y_result_temp(j) = real(reflectivity)^2+ imag(reflectivity)^2;
    j= j+1;
end
%% Plot the strained spectrum
if(usingApp)
    hold(app.FBGStrainedGraph, 'on');
    plot(app.FBGStrainedGraph, currentLarray, y_result_strained);
    plot(app.FBGStrainedGraph, currentLarray, y_result_temp);
    textPosition = [0.95, 0.9]; % Adjust text position as needed
    text(app.FBGStrainedGraph, textPosition(1), textPosition(2), sprintf('Bragg wavelength (temp/temp+strain): %.3f/%.3f nm', newBraggL*1e9, newBraggL_strain*1e9), ...
        'Units', 'Normalized', 'HorizontalAlignment', 'right');
    xlim(app.FBGStrainedGraph, [currentLarray(1) currentLarray(end)]);
    ylim(app.FBGStrainedGraph, "auto");
    hold(app.FBGStrainedGraph, 'off');
    app.simulatedXResp = currentLarray;
    app.simulatedYResp_strain = y_result_strained;
    app.simulatedYResp_temp = y_result_temp;
else
    figure;
    hold on;
    plot(currentLarray, y_result_strained);
    plot(currentLarray, y_result_temp);
end





