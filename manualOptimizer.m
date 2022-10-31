%% Manual Engine optimizer

CDA_f = 0.3*0.00000672*1.5; 
CDA_o = 0.25*0.00002462*1.77;
throatDiam = 1;
A_throat = pi()*(convlength(throatDiam, 'in', 'm')/2)^2; % Throat area in m

% CDA_o = x(1);
% CDA_f = x(2);
% A_throat = x(3);

params = getParams(); 
params.CDA_f = CDA_f;
params.CDA_o = CDA_o;
params.A_throat = A_throat;
params.timestep = 0.012;
totalVolume = pi*((convlength(3.75, 'in', 'm')/2)^2)*(convlength(5, 'ft', 'm')+convlength(10, 'in', 'm'));
OFRatio = 3.8;
params.V_f = (1/(1+OFRatio))*totalVolume; 
params.V_o = (OFRatio/(1+OFRatio))*totalVolume;

[data, avgThrust, maxThrust, impulse, isp] = runEngineSim(params, 1);
disp(sprintf("The total impulse of this engine is %.7f N*s", impulse));
disp(sprintf("The average thrust of this engine is %.7f lbf", avgThrust*0.2248089431));
disp(sprintf("The average chamber pressure of this engine is %.7f psi \n \n", mean(data.P_c)*1.450377e-1));

% Calculate Ideal Expansion Ratio
% Based off Nakka equation: http://www.nakka-rocketry.net/articles/noz_example3.pdf


gam = data.gam_0(50);
term1 = ((gam + 1)/2)^(1/(gam-1)); 
[~, ~, P_e, ~] = atmosisa(params.alt); % Pressure in Pascals
term2 = (P_e/data.P_c(50)/1000)^(1/gam); 
term3 = sqrt(((gam+1)/(gam-1))*(1-(P_e/data.P_c(50)/1000)^((gam-1)/gam))); 
A_ratio = term1 * term2 * term3; % area ratio, throat over exit
optimalExit = throatDiam*(1/sqrt(A_ratio)); % optimal exit diameter, inches






