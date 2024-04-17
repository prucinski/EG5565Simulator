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


    dT = -1;
    strain = 0.0005;

    %good explaination here too: https://www.mdpi.com/1424-8220/11/7/6926

    %% get the lambda term (paper 4, eq 4)
    term0 = (pi*r_f^2/(2*h*r_p*E_h))+1/E_f;
    term1 = 2*r_p/(pi*r_f^2)*term0;
    integralFunction = @(theta)  1./(r_p.*(1-sin(theta))./G_a+r_p./G_p.*log(r_p./r_f));
    term2 = integral(integralFunction, 0, acos(b_rp));
    lambdaTerm = sqrt(term1*term2);
    sectionL = 0.005;
    x = 0;
    x_vals = linspace(0,sectionL, 50);
    y_res = [];
    divisor = cosh(sectionL*lambdaTerm);

    term_e1 = (a_h - a_f)*dT/term0;
    for i= 1:length(x_vals)
        term_e1 = (a_h - a_f)*dT/term0;
        e_t = term_e1*term_e2 + a_f*dT;
        term_e2 = 1 - cosh(lambdaTerm*x_vals(i))/divisor;
        e_m = strain/(E_f*term0)*term_e2;
        y_res(i) = e_m;
    end
    plot(x_vals, y_res);


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
