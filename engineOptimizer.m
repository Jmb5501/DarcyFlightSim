%% Optimizer
% For optimizing for a specified average thrust, will optimize for hitting
% target thrust at the highest total impulse possible

params = getParams();
x0 = [params.CDA_o, params.CDA_f, params.A_throat];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.25*0.00002462*0.1, 0.3*0.00000672*0.1, 6.413016532327705e-04*0.1];
ub = [0.25*0.00002462*8, 0.3*0.00000672*8, 6.413016532327705e-04*8];
% CDA_f = 0.3*0.00000672*2; 
% CDA_o = 0.25*0.00002462*2.4;
% throatDiam = 1.125;
% A_throat = pi()*(convlength(throatDiam, 'in', 'm')/2)^2; % Throat area in m
options = optimoptions("fmincon", "StepTolerance", 1e-08, "Display","iter");

x = fmincon(@func, x0, A, b, Aeq, beq, lb, ub, [], options);

function [val] = func(inputs)
    params = getParams();
    params.CDA_f = inputs(2);
    params.CDA_o = inputs(1);
    params.A_throat = inputs(3); 
    totalVolume = pi*(convlength(3.75, 'in', 'm')^2)*(convlength(5, 'ft', 'm')+convlength(10, 'in', 'm'));
    params.V_f = (1/4.8)*totalVolume; 
    params.V_o = (3.8/4.8)*totalVolume;
    [~, avgThrust, ~, impulse, ~] = runSim(params, 0);
    targetThrust = 2224; % N, target thrust that we are solving for
    val = 10000*abs(targetThrust - avgThrust) - impulse;
end