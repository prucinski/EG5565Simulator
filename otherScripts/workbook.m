        E_h = 72;           %GPa, Young's modulus of host structure
        E_f = 72;           %GPa, Young's modulus of optical fiber
        G_a = 0.714;     %MPa, shear modulus of adhesive
        G_p = 2.02;       %MPa, shear modulus of coating
        r_p = 125;          %micrometers, radius of coating
        r_f = 62.5;         %micrometers, radius of fiber
        h =   0.1;          %ometers, thickness of host material (tank by def 0.1m)

        b_rp = 0.2;         %1, gap between adhesive and fiber 
        a_h = 23;           %microe/K, coeff of thermal expansion of host
        a_f = 0.5;          %microe/K, coeff of thermal expansion of fiber

    %convert to SI units
    G_a = G_a*1e6;   %MPa to Pa     %MPa to GPa
    G_p = G_p*1e6;   %MPa to Pa     %MPa to GPa

    E_h = E_h*1e9;  %GPa to Pa
    E_f = E_f*1e9;  %GPa to Pa

    r_p = r_p*1e-6;  %um to m
    r_f = r_f*1e-6;  %um to m
    a_h = a_h*1e-6; %ue/K to e/K
    a_f = a_f*1e-6; %ue/K to e/K


    dT = 0;
    strain = 0.01;

    %good explaination here too: https://www.mdpi.com/1424-8220/11/7/6926

    %% get the lambda term (paper 4, eq 4)
% Define range of sectionL values
sectionL_values = [0.005, 0.01, 0.025, 0.04, 0.06]; % Modify as needed

% Initialize plot
figure;
hold on;
grid on;

% Loop through sectionL values
for sectionL = sectionL_values
    % Initialize arrays
    x_vals = linspace(-sectionL, sectionL, 100);
    y_res = zeros(size(x_vals));
    
    % Calculate divisor
    divisor = cosh(sectionL * lambdaTerm);
    
    % Calculate strains
    for i = 1:length(x_vals)
        term_e2 = 1 - cosh(lambdaTerm * abs(x_vals(i))) / divisor;
        e_m = strain / (E_f * term0) * term_e2;
        y_res(i) = e_m / strain;
    end
    
    % Plot results
    plot(x_vals, y_res, 'DisplayName', sprintf('bondedL = %.3fm', sectionL));
end

% Add labels and legend
xlabel('Length along the bonded fibre');
ylabel('Normalized Strain');
title('Normalized Strain vs length along fibre for different bonding lenghts');
legend('Location', 'northeast');


    e_t = term_e1*term_e2 + a_f*dT;

    %% get mechanical strain
    e_m = strain/(E_f*term0)*term_e2;
    %disp("Strain transferred " + e_m);


    disp(['Lambda term: ', num2str(lambdaTerm)]);
     disp(['term0: ', num2str(term0)]);
    disp(['cosh(lambdaTerm): ', num2str(cosh(lambdaTerm))]);
    disp(['cosh(x*lambdaTerm): ' , num2str(cosh(x*lambdaTerm))]);
    disp(['cosh(sectionL*lambdaTerm): ' , num2str(cosh(sectionL*lambdaTerm))]);
    disp(['cosh(x*lambda)/cosh(sectionL*lambda)' , num2str(cosh(x*lambdaTerm)/divisor)]);



    %let's try getting the information back (from highest peak)
    recoveredHeight = 1.550063e-06;
    k_e = 0.780096011500;
    e_m_recovered = (recoveredHeight - 1.5500e-06)/(k_e*1.5500e-06);
    disp(e_m_recovered);
    e_actual = e_m_recovered*E_f*term0/(1 - cosh(lambdaTerm*0)/divisor);    %0 is the most powerful component so I hope this can be used
    disp(e_actual);
    %close enough