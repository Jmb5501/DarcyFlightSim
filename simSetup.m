%%%%% Simulation Setup for Darcy II %%%%%
% Purpose: Allows users to specify vehicle parameters, engine model, aero 
% decks, simuation paramters, and initial conditions
%
% THIS IS THE ONLY SCRIPT YOU NEED TO EDIT
% 
% Original Author: John Beavers
% Last Modified: John Beavers (10/20/22)

% TODO - build separate Mass Prop Estimator script and use that to get CGs
% and I

function [vehicle,engine,chute,aero_deck,sim,env,IC] = simSetup()

%% Vehicle Parmeters
% PARAMETERS WITH DISTRIBUTIONS (used for Monte Carlo dispersion analysis)
% Note: nom = nominal, upper = upper bound, lower = lower bound
vehicle.drymass.nom = 11.3398; % wet mass [kg]
    vehicle.wetmass.upper =  110.327;
    vehicle.wetmass.lower = 104.945;

% Wet Centers of Gravity (x,y,z) measured from nosecone tip [m]
vehicle.cgWet_x.nom = -3.5517;
    vehicle.cgWet_x.upper = -3.4;
    vehicle.cgWet_x.lower = -3.6;
vehicle.cgWet_y.nom = 0.0009;
    vehicle.cgWet_y.upper = 0.0010;
    vehicle.cgWet_y.lower = 0.0008;
vehicle.cgWet_z.nom = -0.006;
    vehicle.cgWet_z.upper = -0.005;
    vehicle.cgWet_z.lower = -0.004;

% Dry Centers of Gravity (x,y,z) measured from nosecone tip [m]
vehicle.cgDry_x.nom = -3.482;
    vehicle.cgDry_x.upper = -3.4;
    vehicle.cgDry_x.lower = -3.6;
vehicle.cgDry_y.nom = 0.0011;
    vehicle.cgDry_y.upper = 0.0010;
    vehicle.cgDry_y.lower = 0.0008;
vehicle.cgDry_z.nom = -0.000765;
    vehicle.cgDry_z.upper = -0.005;
    vehicle.cgDry_z.lower = -0.004;

% Wet MOI Tensor [kg*m^2]
vehicle.IxxWet.nom = 0.663228;
    vehicle.IxxWet.upper = 0.673228;
    vehicle.IxxWet.lower = 0.653228;
vehicle.IxyWet.nom = 0.132139;
    vehicle.IxyWet.upper = 0.142139;
    vehicle.IxyWet.lower = 0.122139;
vehicle.IxzWet.nom = 0.06769;
    vehicle.IxzWet.upper = 0.06869;
    vehicle.IxzWet.lower = 0.06669;
vehicle.IyyWet.nom = 181.0307;
    vehicle.IyyWet.upper = 182.0307;
    vehicle.IyyWet.lower = 180.0307;
vehicle.IyzWet.nom = -0.00349;
    vehicle.IyzWet.upper = -0.00359;
    vehicle.IyzWet.lower = -0.00349;
vehicle.IzzWet.nom = 181.0307;
    vehicle.IzzWet.upper = 182.0307;
    vehicle.IzzWet.lower = 180.0307;

% Dry MOI Tensor [kg*m^2]
vehicle.IxxDry.nom = 0.611089;
    vehicle.IxxDry.upper = 0.623228;
    vehicle.IxxDry.lower = 0.603228;
vehicle.IxyDry.nom = 0.139071;
    vehicle.IxyDry.upper = 0.142139;
    vehicle.IxyDry.lower = 0.122139;
vehicle.IxzDry.nom = 0.072565;
    vehicle.IxzDry.upper = 0.07469;
    vehicle.IxzDry.lower = 0.06669;
vehicle.IyyDry.nom = 173.2412;
    vehicle.IyyDry.upper = 173.2412;
    vehicle.IyyDry.lower = 173.2412;
vehicle.IyzDry.nom = -0.0035;
    vehicle.IyzDry.upper = -0.00359;
    vehicle.IyzDry.lower = -0.00349;
vehicle.IzzDry.nom = 173.2412;
    vehicle.IzzDry.upper = 182.0307;
    vehicle.IzzDry.lower = 180.0307;

% PARAMETRS WITH KNOWN CONSTANT VALUES   
vehicle.diameter = 0.1016; % diameter [m]
vehicle.length = 3.048;% length [m]

%% Engine Model
% runs the engine simulator to generate thrust curve and propellant masses
% cant - engine cant angle [deg] (1st index is angle w.r.t body x axis, 2nd index is w.r.t y axis)

engine_params = engineSetup(); 

% Code currently copied from manualOptimizer.m
engine_params.CDA_f =  0.3*0.00000672*1.5; 
engine_params.CDA_o = 0.25*0.00002462*1.77;

[engine, ~, ~, ~, ~] = runEngineSim(engine_params, false);

%% Parachute Model
% main_mass - mass of main chute [kg]
% main_CD - drag coefficient of main chute
% main_A - area of main chute [m^2]
% main_deploy - height above initial launch altitude to deploy main chute [m] 
% drogue_mass - mass of drogue chute [kg]
% drogue_CD - drag coefficient of drogue chute
% drogue_A - area of drogue chute [m^2]
% DROGUE CHUTE IS HARDCODED IN TO DEPLOY AT APOGEE

chute.main_mass = 0.380;
chute.main_CD = 2.2;
chute.main_A = 2.6263689; % 72" chute
chute.main_deploy = 1000 / 3.281;
chute.drogue_mass = 0.062;
chute.drogue_CD = 1.5;
chute.drogue_A = 0.29171555; % 24" elliptical drogue

%% Aero Deck 
% specify the aero deck to be used
aero_deck = readtable('Darcy2_Aero.csv',"FileType",'spreadsheet'); 

%% Simulation Parameters
% dofs - specify what DOF mode to run the sim in (1,3,6)
% trials - number of trials to run
% method - 'mc' for Monte Carlo, 'sweep_{parameter}' to sweep over a single vehicle parameter, 'single' for a single simulation
% dt - simulation time step [s]
% tspan - simulation time [s]
% flag - set to 1 to run until apogee is hit, 2 to run until end of sim time or impact
% integrator - 1 for RK4, 2 for Euler

sim.dofs = 1;
sim.trials = 1;
sim.method = 'single';
sim.dt = .1; % MAKE SURE THIS dt AND THE ENGINE SIM dt ARE EQUAL (check getParams.m for engine dt)
sim.tspan = 1000;
sim.flag = 2;
sim.integrator= 1;

if sim.dt ~= engine_params.timestep
    error('Flight Sim dt and Engine Sim dt must be equal! Exiting program...')
end
%% Environmental Parameters
% gamma - ratio of specific heat for air
% R - gas constant for air [J/(kg*K)]
% Re - radius of Earth [m]
% g0 - gravitational acceleration at sea level (m/s^2)


env.gamma = 1.4;
env.R = 287;
env.Re = 6.3781E6; 
env.g= 9.81;
%% Initial Condtions
% LLA - latitude [deg], longitude [deg], altitude [m]
% pos - position in NED frame [m]
% vel - velocity in NED frame [m/s]
% roll - initial roll angle of vehicle when on launch rail [deg]
% rail_angle - angle of launch rail w.r.t the D axis of the NED frame [deg]
%              1st index is elevation (measured from zenith, positive rotation along
%              negative E axis), 2nd index is azimuth (measured cw from north, positive rotation along positive D axis)
% rail_length - length of launch rail [m]


IC.LLA = [30.933678; -81.518268; 3]; 
IC.pos = [0;0;0];
IC.vel = [0;0;0];
IC.roll = 0;
IC.rail_angle = [.5,0];
IC.rail_length = 10;
end