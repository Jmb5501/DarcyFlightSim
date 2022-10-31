%%%%% 3DOF Trajectory Model for Darcy II%%%%%
% Purpose: 3DOF (3D position) simulator for YJSP-XL vehicles
% Key Features:
%   1) Model includes gravity, thrust, reco, and aerodynamic force in axial direction.
%   2) Integrator is built-in RK4
%   3) Capable of modeling effects due to launch rail angle and wind.
% Original Author: John Beavers
% Last Modified: John Beavers
% Date Updated: 10/24/22

function flight_data = traj_3dof(Scenario,aero_deck)

%% Unpack Scenario Structure (only what is needed for 3DOF simulation)
drymass = Scenario.drymass;
diameter = Scenario.diameter;
mdot = Scenario.mdot; % vector
thrust = Scenario.thrust; % vector
propmass = Scenario.propmass;
main_mass = Scenario.main_mass;
main_CD = Scenario.main_CD;
main_A = Scenario.main_A;
main_deploy = Scenario.main_deploy;
drogue_mass = Scenario.drogue_mass;
drogue_CD = Scenario.drogue_CD;
drogue_A = Scenario.drogue_A;
env = Scenario.env;
pos_0 = Scenario.pos;
vel_0 = Scenario.vel;
alt_0 = Scenario.LLA(3);
rail_angle = Scenario.rail_angle;
dt = Scenario.dt;
tf = Scenario.tspan;
flag = Scenario.flag;
integrator = Scenario.integrator;
%% Time Definitions
apogee_flag = 0;
main_deploy_flag = 0;
meco_flag = 0;
tspan = 0:dt:tf;

%% Propulsion Definitions
mdot_curve = zeros(length(tspan),1);
thrust_curve = zeros(length(tspan),1);
mdot_curve(1:length(mdot)) = mdot;
thrust_curve(1:length(thrust)) = thrust;

%% Initialization
x = [pos_0',vel_0',drymass+propmass+main_mass+drogue_mass]; 
flight_data = zeros(length(tspan),length(x)+3); % three other 0s are for Mach, Drag, and Thrust outputting

%% Integration (RK4 or Euler)
for i = 1:length(tspan)
    if integrator == 1
        x_next = RK4(dt,x,diameter,mdot_curve(i),thrust_curve(i),main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    elseif integrator == 2
        x_next = euler(dt,x,diameter,mdot_curve(i),thrust_curve(i),main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    else
        error('Invalid input for integrator, check sim_setup! Exiting program...')
    end
    flight_data(i,1:7) = x_next;
    x = x_next;
    flight_data(i,[3,6]) = -flight_data(i,[3,6]); % since we are in NED frame, negative positions are "up", this corrects this for visualization and logic purposes
    
    % Copied code from rockdetdyn_3dof to output Mach, Drag, and Thrust
    vel_mag = norm([x_next(4),x_next(5),x_next(6)]);
    [rho,a] = stdatmo(-x_next(3)+alt_0);
    mach = vel_mag/a;
    alpha = 0;
    if mach == 0
        CD = 0;
    else
        CD = interp1(aero_deck{aero_deck{:,2}==alpha,1},aero_deck{aero_deck{:,2}==alpha,3},mach);
    end
    area = pi*(diameter/2)^2;
    D_aero = 1/2*rho*vel_mag^2*area*CD;
    flight_data(i,8:10) = [mach,D_aero,thrust_curve(i)];

    % logic for displaying when MECO is hit
    if thrust_curve(i) == 0 && ~meco_flag
        meco_flag = 1;
        disp(['MECO at ' num2str(flight_data(i,3)*3.2804) ' ft']);
    end
    
    % logic for displaying when drogue chute has been deployed in console
    if i>1 && flight_data(i,6) <= 0 && ~apogee_flag
        apogee_flag = 1;
        disp(['Drogue chute delpoyed at ' num2str(flight_data(i,3)*3.2804) ' ft']);
    end
    
    % logic for displaying when main chute has been deployed in console
    if flight_data(i,3) <= main_deploy && apogee_flag && ~main_deploy_flag
        main_deploy_flag = 1;
        disp(['Main chute delpoyed at ' num2str(flight_data(i,3)*3.2804) ' ft']);
        disp(['Drogue chute descent rate:  ' num2str(flight_data(i-1,6)*3.2804) ' ft/s']);
    end

    
    % logic for checking if rocket has hit ground
    if i>1 && flight_data(i,3) <= 0
        break
    end    

    % logic for ending sim when apogee is found
    if flag == 1  && apogee_flag
        break   
    else
        continue
    end
end


end

%%% HELPER FUNCTIONS %%%

%% RK4 Integrator
% Purpose: RK4 integrator for simulation of a rocket (body frame)
function [x_next] = RK4(dt,x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck)
    y1 = x;
    k1 = rocketdyn_3dof(y1,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);

    y2 = x + 1/2*k1*dt;
    k2 = rocketdyn_3dof(y2,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    
    y3 = x + 1/2*k2*dt;
    k3 = rocketdyn_3dof(y3,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    
    y4 = x + k3*dt;
    k4 = rocketdyn_3dof(y4,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    
    x_next = x + 1/6*(k1 + 2*k2 + 2*k3 + k4)*dt;
end
%% Euler Integrator
% Purpose: Euler integrator for simulation of a rocket
function [x_next] = euler(dt,x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck)

    dxdt = rocketdyn_3dof(x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck);
    
    x_next = x + dxdt*dt;

end
%% 3-DOF Rocket Dynamics Model %%
% Purpose: State space physics model for 3-DOF simulations of a rocket
function [dxdt] = rocketdyn_3dof(x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,rail_angle,apogee_flag,aero_deck)
    
    % States - NED frame
    pos_x = x(1);
    pos_y = x(2);
    pos_z = x(3);
    vel_x = x(4);
    vel_y = x(5);
    vel_z = x(6);
    mass = x(7);
    
    % Launch Rail DCM - Body to NED frame DCM
    el = rail_angle(1);
    az = rail_angle(2);
    el_rot = [cosd(el), 0, -sind(el);
                 0     ,1,    0;
              sind(el), 0, cosd(el)];

    az_rot = [cosd(az), -sind(az), 0;
              sind(az),  cosd(az), 0;
                0       ,   0,     1];

    vec = [0; 0; -1]; %initial launch rail angle is pointed vertically
    launch_vec = az_rot*el_rot*vec; %rotate to get launch vector

    % T = engine thrust model
%     engine_cant_DCM_Z = [1,0,0;
%                         0, cosd(engine_cant(1)), -sind(engine_cant(1));
%                         0, sind(engine_cant(1)), cosd(engine_cant(1))]; % rotation about body z axis
%     engine_cant_DCM_X = [cosd(engine_cant(2)), -sind(engine_cant(2)), 0;
%                          sind(engine_cant(2)), cosd(engine_cant(2)), 0;
%                          0, 0, 1]; % rotation about body x axis
%     thrust_unit = engine_cant_DCM_X*engine_cant_DCM_Z*[-1;0;0];
    T = thrust*launch_vec; % thrust in body frame

    % D_aero = aerodynamic drag force model
    vel_mag = norm([vel_x,vel_y,vel_z]);
    if vel_mag == 0
        vel_unit = [0,0,0];
    else
        vel_unit = [vel_x,vel_y,vel_z]/vel_mag;
    end
    area = pi*(diameter/2)^2;
    [rho,a] = stdatmo(-pos_z+alt_0); % compute air density and speed of sound at current altitude
    mach = vel_mag/a; % mach number
    alpha = 0; % angle of attack (deg)
    if mach == 0
        CD = 0;
    else
        % linear interpolation of Mach number to get CD
        CD = interp1(aero_deck{aero_deck{:,2}==alpha,1},aero_deck{aero_deck{:,2}==alpha,3},mach);
    end
    D_aero = (1/2)*rho*vel_mag^2*area*CD*-vel_unit;

    %D_reco = reco drag force model
    if apogee_flag
        D_reco = (1/2)*rho*vel_mag^2*drogue_A*drogue_CD*-vel_unit;
        if -pos_z <= main_deploy
            D_reco = (1/2)*rho*vel_mag^2*(drogue_A*drogue_CD+main_A*main_CD)*-vel_unit;
        end
    else
        D_reco = [0,0,0];
    end

    % W = rocket weight model
    if pos_z == 0
        W = [0,0,0];
    else
        W = [0,0,env.g*mass*(env.Re/(-pos_z+env.Re))^2];
    end

    
    dxdt = zeros(size(x));
    dxdt(1) = vel_x;
    dxdt(2) = vel_y;
    dxdt(3) = vel_z;
    dxdt(4) = (T(1)+D_aero(1)+D_reco(1)+W(1))/mass;
    dxdt(5) = (T(2)+D_aero(2)+D_reco(2)+W(2))/mass;
    dxdt(6) = (T(3)+D_aero(3)+D_reco(3)+W(3))/mass;
    dxdt(7) = -mdot;
end