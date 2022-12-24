function stabilityPlot(flight_data,dofs,dt,cg_wet,cg_dry,aero_deck,diameter)
%%%%% Stability Plotter for Darcy II %%%%%
% Purpose: Plots stability profile for 1DOF, 3DOF, and 6DOF Darcy II sims
% Original Author: John Beavers
% Last Modified: John Beavers
% Date Updated: 12/06/22

% Intialization
t = 0:dt:(length(flight_data)-1)*dt;
cg_data = (-0.0019.*t.^3- 0.0182.*t.^2 +0.7622.*t+ 87.27)/39.37;
mask = t >=16.562;
cg_data(mask) = cg_data(t == 16.562);
cp_data = zeros(length(flight_data),1);
static_margin = zeros(length(flight_data),1);
mass_data = flight_data(:,3);
alpha = 0;

% Loop to calcuate stability data
for i = 1:length(mass_data)
    mach = flight_data(i,4);
    cp_data(i) = (interp1(aero_deck{aero_deck{:,2}==alpha,1},aero_deck{aero_deck{:,2}==alpha,13},mach))/39.37; %inches to m conversion
    
    if mach < 0.1
        static_margin(i) = NaN;
    else
        static_margin(i) = (cp_data(i) - cg_data(i))/diameter;
    end
end

% Plotting
figure
plot(t,static_margin)
title("Stability Plot")
ylabel("Cal")
xlabel("Time (s)")