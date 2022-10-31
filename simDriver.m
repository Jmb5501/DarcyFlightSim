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

        % Post Processing
        mask = [true;flight_data(2:end,1) ~= 0];
        flight_data = flight_data(mask,:);
        disp(['Main chute descent rate:  ' num2str(flight_data(end,2)*3.2804) ' ft/s']);
    
        % Plotting
        time = 0:sim.dt:(length(flight_data)-1)*sim.dt;
        subplot(3,3,1)
        plot(time,flight_data(:,1)*3.281) % convert m to ft
        grid
        xlabel('Time (s)')
        ylabel('Position (ft)')
        subplot(3,3,2)
        plot(time,flight_data(:,2)*3.281) % convert m/s to ft/s
        grid
        xlabel('Time (s)')
        ylabel('Velocity (ft/s)')
        subplot(3,3,3)
        plot(time,flight_data(:,4))
        grid
        xlabel('Time (s)')
        ylabel('Mach')
        subplot(3,3,4)
        plot(time,flight_data(:,5)/4.448) % convert N to lbf
        grid
        xlabel('Time (s)')
        ylabel('Aerodynamic Drag (lbf)')
        subplot(3,3,5)
        plot(time,flight_data(:,6)/4.448) % convert N to lbf
        grid
        xlabel('Time (s)')
        ylabel('Thrust (lbf)')
        subplot(3,3,6)
        plot(time,flight_data(:,3)*2.205) % convert kg to lb
        grid
        xlabel('Time (s)')
        ylabel('Mass (lb)')
        subplot(3,3,7)
        plot(time(1:end-1),(diff(flight_data(:,2))./diff(time))/9.81) % convert m/s^2 to g's
        grid
        xlabel('Time (s)')
        ylabel('Acceleration (g)')

    elseif sim.dofs == 3
        flight_data = traj_3dof(Scenarios,aero_deck);

        % Post Processing
        mask = any([true*ones(1,3);flight_data(2:end,4:6) ~= 0],2);
        flight_data = flight_data(mask,:);
        disp(['Main chute descent rate:  ' num2str(flight_data(end,6)*3.2804) ' ft/s']);
    
        % Plotting
        time = 0:sim.dt:(length(flight_data(:,1))-1)*sim.dt;
        subplot(3,3,1)
        plot(time,flight_data(:,3)*3.281) % convert m to ft
        grid
        xlabel('Time (s)')
        ylabel('Position (ft)')
        subplot(3,3,2)
        vel_mags = (flight_data(:,4).^2+flight_data(:,5).^2+flight_data(:,6).^2).^(1/2);
        plot(time,vel_mags*3.281) % convert m/s to ft/s
        grid
        xlabel('Time (s)')
        ylabel('Velocity (ft/s)')
        subplot(3,3,3)
        plot(time,flight_data(:,8))
        grid
        xlabel('Time (s)')
        ylabel('Mach')
        subplot(3,3,4)
        plot(time,flight_data(:,9)/4.448) % convert N to lbf
        grid
        xlabel('Time (s)')
        ylabel('Aerodynamic Drag (lbf)')
        subplot(3,3,5)
        plot(time,flight_data(:,10)/4.448) % convert N to lbf
        grid
        xlabel('Time (s)')
        ylabel('Thrust (lbf)')
        subplot(3,3,6)
        plot(time,flight_data(:,7)*2.205) % convert kg to lb
        grid
        xlabel('Time (s)')
        ylabel('Mass (lb)')
        subplot(3,3,7)
        plot(time(1:end-1),(diff(vel_mags)./diff(time))/9.81) % convert m/s^2 to g's
        grid
        xlabel('Time (s)')
        ylabel('Acceleration (g)')
        
        % trajectory plot
        figure
        plot3(flight_data(:,1)*3.281,flight_data(:,2)*3.281,flight_data(:,3)*3.281) % convert m to ft
        grid
        xlabel('X Position(ft)')
        ylabel('Y Position(ft)')
        zlabel('Z Position(ft)')

    elseif sim.dofs == 6
        flight_data = traj_6dof(Scenarios,aero_deck);
    else
        error('Input for sim.dofs is invalid, check sim_setup for a valid DOF entry! Exiting program...')
    end
    

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
