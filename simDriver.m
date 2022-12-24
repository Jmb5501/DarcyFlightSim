%%%%% Trajectory Simulator Driver Script %%%%%
% Purpose: Driver Script for Executing the YJSP_XL Trajectory Simulator
% THIS IS THE ONLY SCRIPT YOU NEED TO RUN!
% Key Notes:
%   1) Can run 1 DOF (vertical position only), 3 DOF (3D position only) or 6 DOF sims
%   2) Can run Monte Carlo (MC) simulations
% 
% Original Author: John Beavers
% Last Modified: John Beavers (10/24/22)


clear
clc
close all
tic
% Import data from sim_setup script
[vehicle,engine,chute,aero_deck,sim,env,IC] = simSetup();

% Initialize Scenarios structure based on how many sims to run
Scenarios(sim.trials) = struct();

%% Single Simulation
if strcmp(sim.method,'single')
    Scenarios.drymass = vehicle.drymass.nom;
    Scenarios.cgWet_x = vehicle.cgWet_x.nom;
    Scenarios.cgWet_y = vehicle.cgWet_y.nom;
    Scenarios.cgWet_z = vehicle.cgWet_z.nom;
    Scenarios.cgDry_x = vehicle.cgDry_x.nom;
    Scenarios.cgDry_y = vehicle.cgDry_y.nom;
    Scenarios.cgDry_z = vehicle.cgDry_z.nom;
    Scenarios.IxxWet = vehicle.IxxWet.nom;
    Scenarios.IxyWet = vehicle.IxyWet.nom;
    Scenarios.IyyWet = vehicle.IyyWet.nom;
    Scenarios.IyzWet = vehicle.IyzWet.nom;
    Scenarios.IzzWet = vehicle.IzzWet.nom;
    Scenarios.IxxDry = vehicle.IxxDry.nom;
    Scenarios.IxyDry = vehicle.IxyDry.nom;
    Scenarios.IyyDry = vehicle.IyyDry.nom;
    Scenarios.IyzDry = vehicle.IyzDry.nom;
    Scenarios.IzzDry = vehicle.IzzDry.nom;
    Scenarios.diameter = vehicle.diameter;
    Scenarios.length = vehicle.length;
    Scenarios.mdot = engine.m_dot;
    Scenarios.thrust = engine.thrust;
    Scenarios.propmass = engine.m_f(1)+engine.m_N2O(1);
    Scenarios.main_mass = chute.main_mass;
    Scenarios.main_CD = chute.main_CD;
    Scenarios.main_A = chute.main_A;
    Scenarios.main_deploy = chute.main_deploy;
    Scenarios.drogue_mass = chute.drogue_mass;
    Scenarios.drogue_CD = chute.drogue_CD;
    Scenarios.drogue_A = chute.drogue_A;
    Scenarios.dt = sim.dt;
    Scenarios.tspan = sim.tspan;
    Scenarios.flag = sim.flag;
    Scenarios.integrator = sim.integrator;
    Scenarios.env = env;
    Scenarios.LLA = IC.LLA;
    Scenarios.pos = IC.pos;
    Scenarios.vel = IC.vel;
    Scenarios.roll = IC.roll;
    Scenarios.rail_angle = IC.rail_angle;
    Scenarios.rail_length = IC.rail_length;
    
    if sim.dofs == 1
        flight_data = traj_1dof(Scenarios,aero_deck);
    elseif sim.dofs == 3
        flight_data = traj_3dof(Scenarios,aero_deck);
    elseif sim.dofs == 6
        flight_data = traj_6dof(Scenarios,aero_deck);
    else
        error('Input for sim.dofs is invalid, check sim_setup for a valid DOF entry! Exiting program...')
    end
    flightProfilePlot(flight_data,sim.dofs,sim.dt)
    stabilityPlot(flight_data,sim.dofs,sim.dt,Scenarios.cgWet_x,Scenarios.cgDry_x,aero_deck,Scenarios.diameter)

%% Monte Carlo
elseif strcmp(sim.method,'mc')
    % TODO: Add in Monte Carlo math for scenario creation
    for i = 1:sim.trials
        Scenarios(i).wetmass = vehicle.wetmass.nom;
        Scenarios(i).propmass = vehicle.propmass.nom;
        Scenarios(i).cgWet_x = vehicle.cgWet_x.nom;
        Scenarios(i).cgWet_y = vehicle.cgWet_y.nom;
        Scenarios(i).cgWet_z = vehicle.cgWet_z.nom;
        Scenarios(i).cgDry_x = vehicle.cgDry_x.nom;
        Scenarios(i).cgDry_y = vehicle.cgDry_y.nom;
        Scenarios(i).cgDry_z = vehicle.cgDry_z.nom;
        Scenarios(i).IxxWet = vehicle.IxxWet.nom;
        Scenarios(i).IxyWet = vehicle.IxyWet.nom;
        Scenarios(i).IyyWet = vehicle.IyyWet.nom;
        Scenarios(i).IyzWet = vehicle.IyzWet.nom;
        Scenarios(i).IzzWet = vehicle.IzzWet.nom;
        Scenarios(i).IxxDry = vehicle.IxxDry.nom;
        Scenarios(i).IxyDry = vehicle.IxyDry.nom;
        Scenarios(i).IyyDry = vehicle.IyyDry.nom;
        Scenarios(i).IyzDry = vehicle.IyzDry.nom;
        Scenarios(i).IzzDry = vehicle.IzzDry.nom;
        Scenarios(i).diameter = vehicle.diameter;
        Scenarios(i).length = vehicle.length;
        Scenarios(i).dt = sim.dt;
        Scenarios(i).tspan = sim.tspan;
        Scenarios(i).gamma = env.gamma;
        Scenarios(i).R = env.R;
        Scenarios(i).g = env.g;
        Scenarios(i).rail_angle = env.rail_angle;
        Scenarios(i).rail_length = env.rail_length;
        Scenarios(i).LLA = IC.LLA;
        Scenarios(i).pos = IC.pos;
        Scenarios(i).vel = IC.vel;
        Scenarios(i).roll = IC.roll;
        
        if sim.dofs == 1
            flight_data = traj_1dof(Scenarios(i),aero_deck);

            % Post Processing

        elseif sim.dofs == 3
            flight_data = traj_3dof(Scenarios(i),aero_deck);
        elseif sim.dofs == 6
            flight_data = traj_6dof(Scenarios(i),aero_deck);
        else
            error('Input for sim.dofs is invalid, check sim_setup for a valid DOF entry! Exiting program...')
        end
    end



%% Single Variable Sweep
elseif contains(sim.method,'sweep')



%% Error Handle
else
    error('Input for sim.method is not valid, check sim_setup for a valid method entry! Exiting program...')
end
toc