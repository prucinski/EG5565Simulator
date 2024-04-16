%h - host structure, f - fiber, p - protective coating, a - adhesive
% I have no clue what the units for sectionL, fiberL are meant to be so I
% will assume m
function [e_t, e_m] = getThermalAndMechanical(dT, strain, fiberL, sectionL,app) 

    L_f = 0.5e6*fiberL;   %convert to micrometers, half the bonding length (length of fibre)
    S_f = sectionL*1e6;   %convert to micrometers, section length

    %this is a nasty-looking block for detection of errors in the table
    %tried my best to find a programatic solution but to no avail
    %might chuck it into a separate function eventually
    try
      E_f = str2double(app.UITable.Data{1,2});
        if(isnan(E_f)) app.UITable.Data{1,1} = {'72'}; error('E_f is not a number');   end
      E_h = str2double(app.UITable.Data{1,4});
        if(isnan(E_h)) app.UITable.Data{1,4} = {'72'}; error('E_h is not a number');   end
     
     G_p = str2double(app.UITable.Data{3,3});
        if(isnan(G_p)) app.UITable.Data{3,3} = {'2.02'}; error('G_p is not a number');   end
     G_a = str2double(app.UITable.Data{3,4});
        if(isnan(G_a)) app.UITable.Data{3,4} = {'0.714'}; error('G_a is not a number');   end
     
     r_f = str2double(app.UITable.Data{4,2});
        if(isnan(r_f)) app.UITable.Data{4,2} = {'62.5'}; error('r_f is not a number');   end     
     r_p = str2double(app.UITable.Data{4,3});
        if(isnan(r_p)) app.UITable.Data{4,3} = {'125'}; error('r_p is not a number');   end
     b_rp = str2double(app.UITable.Data{4,4});
        if(isnan(b_rp)) app.UITable.Data{4,4} = {'0.2'}; error('b_rp is not a number');   end
     h = str2double(app.UITable.Data{4,5});
        if(isnan(h)) app.UITable.Data{4,5} = {'0.1'}; error('h is not a number');   end
     
     a_f = str2double(app.UITable.Data{5,2});
        if(isnan(a_f)) app.UITable.Data{5,2} = {'0.5'}; error('a_f is not a number');   end
     a_h = str2double(app.UITable.Data{5,5});
        if(isnan(a_h)) app.UITable.Data{5,5} = {'23'}; error('a_h is not a number');   end


    catch ME 
        errordlg("Problem with a variable in the advanced settings tab:  " + ...
                    ME.message + ".Assuming default for variable ", "Data error");
        E_h = 72;           %GPa, Young's modulus of host structure
        E_f = 72;           %GPa, Young's modulus of optical fiber
        G_a = 0.714;     %MPa, shear modulus of adhesive
        G_p = 2.2;       %MPa, shear modulus of coating
        r_p = 125;          %micrometers, radius of coating
        r_f = 62.5;         %micrometers, radius of fiber
        h =   0.1;          %ometers, thickness of host material (tank by def 0.1m)

        b_rp = 0.2;         %1, gap between adhesive and fiber 

        a_h = 23;           %microe/K, coeff of thermal expansion of host
        a_f = 0.5;          %microe/K, coeff of thermal expansion of fiber
    end
    %convert to consistent units
    h = h*1e6;        %m to um
    G_a = G_a*1e-3;   %MPa to GPa
    G_p = G_p*1e-3;   %MPa to GPa
        
    %good explaination here too

    %% get the lambda term (paper 4, eq 4)
    term0 = (pi*r_f^2/(2*h*r_p*E_h))+1/E_f;
    term1 = 2*r_p/(pi*r_f^2)*term0;
    integralFunction = @(theta)  1./(r_p.*(1-sin(theta))./G_a+r_p./G_p.*log(r_p./r_f));
    term2 = integral(integralFunction, 0, acos(b_rp));
    lambdaTerm = sqrt(term1*term2);
    disp(lambdaTerm);
    %% get thermal strain
    %This paper too: https://sensors.myu-group.co.jp/sm_pdf/SM1255.pdf
    %God knows the units they use... Nobody cares about reproducibility
    %anymore
    term_e1 = (a_h - a_f)*dT/term0;
    term_e2 = (1 - (cosh(lambdaTerm*S_f))/(cosh(lambdaTerm*L_f)));
    e_t = term_e1*term_e2 + a_f*dT;

    %% get mechanical strain
    e_m = strain/term0*term_e2;
    disp("Strain transferred " + e_m);
end
