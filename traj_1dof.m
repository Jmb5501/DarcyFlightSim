%%%%% 1DOF Trajectory Model for Darcy II%%%%%
% Purpose: 1DOF (vertical position) simulator for YJSP-XL vehicles
% Key Features:
%   1) Model includes gravity, thrust, reco, and aerodynamic force in axial direction.
%   2) Integrator is built-in RK4
% Original Author: John Beavers
% Last Modified: John Beavers
% Date Updated: 10/24/22

function flight_data = traj_1dof(Scenario,aero_deck)

%% Unpack Scenario Structure (only what is needed for 1DOF simulation)
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
pos_0 = Scenario.pos(3);
vel_0 = Scenario.vel(3);
alt_0 = Scenario.LLA(3);
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
x = [pos_0,vel_0,drymass+propmass+main_mass+drogue_mass]; 
flight_data = zeros(length(tspan),length(x)+3); % three other 0s are for Mach, Drag, and Thrust outputting

%% Integration (RK4 or Euler)
for i = 1:length(tspan)
    if integrator == 1
        x_next = RK4(dt,x,diameter,mdot_curve(i),thrust_curve(i),main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    elseif integrator == 2
        x_next = euler(dt,x,diameter,mdot_curve(i),thrust_curve(i),main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    else
        error('Invalid input for integrator, check sim_setup! Exiting program...')
    end
    flight_data(i,1:3) = x_next;
    x = x_next;
    
    % Copied code from rockdetdyn_1dof to output Mach, Drag, and Thrust
    [rho,a] = stdatmo(x_next(1)+alt_0);
    mach = abs(x_next(2))/a;
    if mach == 0
        CD = 0;
    else
        CD = interp1(aero_deck{aero_deck{:,2}==0,1},aero_deck{aero_deck{:,2}==0,3},mach);
    end
    area = pi*(diameter/2)^2;
    D_aero = 1/2*rho*abs(x_next(2))^2*area*CD;
    flight_data(i,4:6) = [mach,D_aero,thrust_curve(i)];

    % logic for displaying when MECO is hit
    if thrust_curve(i) == 0 && ~meco_flag
        meco_flag = 1;
        disp(['MECO at ' num2str(flight_data(i,1)*3.2804) ' ft']);
    end
    
    % logic for displaying when drogue chute has been deployed in console
    if i>1 && flight_data(i,2) <= 0 && ~apogee_flag
        apogee_flag = 1;
        disp(['Drogue chute delpoyed at ' num2str(flight_data(i,1)*3.2804) ' ft']);
    end
    
    % logic for displaying when main chute has been deployed in console
    if flight_data(i,1) <= main_deploy && apogee_flag && ~main_deploy_flag
        main_deploy_flag = 1;
        disp(['Main chute delpoyed at ' num2str(flight_data(i,1)*3.2804) ' ft']);
        disp(['Drogue chute descent rate:  ' num2str(flight_data(i-1,2)*3.2804) ' ft/s']);
    end

    
    % logic for checking if rocket has hit ground
    if i>1 && flight_data(i,1) <= 0
        break
    end    

    % logic for ending sim when apogee is found
    if flag ==1  && apogee_flag
        break   
    else
        continue
    end
end


end

%%% HELPER FUNCTIONS %%%

%% RK4 Integrator
% Purpose: RK4 integrator for simulation of a rocket
function [x_next] = RK4(dt,x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck)
    y1 = x;
    k1 = rocketdyn_1dof(y1,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);

    y2 = x + 1/2*k1*dt;
    k2 = rocketdyn_1dof(y2,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    
    y3 = x + 1/2*k2*dt;
    k3 = rocketdyn_1dof(y3,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    
    y4 = x + k3*dt;
    k4 = rocketdyn_1dof(y4,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    
    x_next = x + 1/6*(k1 + 2*k2 + 2*k3 + k4)*dt;
end
%% Euler Integrator
% Purpose: Euler integrator for simulation of a rocket
function [x_next] = euler(dt,x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck)

    dxdt = rocketdyn_1dof(x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck);
    
    x_next = x + dxdt*dt;

end
%% 1-DOF Rocket Dynamics Model %%
% Purpose: State space physics model for 1-DOF RK4 simulations of a rocket
function [dxdt] = rocketdyn_1dof(x,diameter,mdot,thrust,main_CD,main_A,drogue_CD,drogue_A,main_deploy,env,alt_0,apogee_flag,aero_deck)
    pos = x(1);
    vel = x(2);
    mass = x(3);

    % T = engine thrust model
    T = thrust;

    % D_aero = aerodynamic drag force model
    area = pi*(diameter/2)^2;
    [rho,a] = stdatmo(pos+alt_0); % compute air density and speed of sound at current altitude
    mach = abs(vel)/a; % mach number
    alpha = 0; % angle of attack (deg)
    if mach == 0
        CD = 0;
    else
        % linear interpolation of Mach number to get CD
        CD = interp1(aero_deck{aero_deck{:,2}==alpha,1},aero_deck{aero_deck{:,2}==alpha,3},mach);
    end
    D_aero = -sign(vel)*(1/2)*rho*abs(vel)^2*area*CD;

    %D_reco = reco drag force model
    if apogee_flag
        D_reco = -sign(vel)*(1/2)*rho*abs(vel)^2*drogue_A*drogue_CD;
        if pos <= main_deploy
            D_reco = -sign(vel)*(1/2)*rho*abs(vel)^2*(drogue_A*drogue_CD+main_A*main_CD);
        end
    else
        D_reco = 0;
    end

    % W = rocket weight model
    if pos == 0
        W = 0;
    else
        W = -env.g*mass*(env.Re/(pos+env.Re))^2;
    end
    
    dxdt = zeros(size(x));
    dxdt(1) = vel;
    dxdt(2) = (T+D_aero+D_reco+W)/mass;
    dxdt(3) = -mdot;
end