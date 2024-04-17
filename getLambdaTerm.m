%h - host structure, f - fiber, p - protective coating, a - adhesive
% I have no clue what the units for sectionL, fiberL are meant to be so I
% will assume m
function [term0, lambdaTerm] = getLambdaTerm(E_f, E_h, G_p, G_a, r_f, r_p, b_rp, h) 

    %convert to SI units
    G_a = G_a*1e6;   %MPa to Pa     
    G_p = G_p*1e6;   %MPa to Pa     
    E_h = E_h*1e9;  %GPa to Pa
    E_f = E_f*1e9;  %GPa to Pa

    r_p = r_p*1e-6;  %um to m
    r_f = r_f*1e-6;  %um to m
        
    %good explaination here too: https://www.mdpi.com/1424-8220/11/7/6926
    disp(E_h);
    %% get the lambda term (paper 4, eq 4)
    term0 = (pi*r_f^2/(2*h*r_p*E_h))+1/E_f;
    term1 = 2*r_p/(pi*r_f^2)*term0;
    integralFunction = @(theta)  1./(r_p.*(1-sin(theta))./G_a+r_p./G_p.*log(r_p./r_f));
    term2 = integral(integralFunction, 0, acos(b_rp));
    lambdaTerm = sqrt(term1*term2);
    disp(['Lambda term: ' , num2str(lambdaTerm)]);
    disp(['term0: ', num2str(term0)]);