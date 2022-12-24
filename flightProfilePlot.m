function flightProfilePlot(flight_data,dofs,dt)
%%%%% Flight Profile Plotter for Darcy II %%%%%
% Purpose: Plots flight profile for 1DOF, 3DOF, and 6DOF Darcy II sims
% Original Author: John Beavers
% Last Modified: John Beavers
% Date Updated: 11/03/22


if dofs == 1
        time = 0:dt:(length(flight_data)-1)*dt;
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
elseif dofs == 3
        time = 0:dt:(length(flight_data(:,1))-1)*dt;
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
        
        % 3D trajectory plot
        figure
        plot3(flight_data(:,1)*3.281,flight_data(:,2)*3.281,flight_data(:,3)*3.281) % convert m to ft
        grid
        xlabel('X Position(ft)')
        ylabel('Y Position(ft)')
        zlabel('Z Position(ft)')
elseif dofs == 6

else
     error('Input for dofs is invalid, unable to plot flight profile! Check sim_setup for a valid DOF entry! Exiting program...')
end
end