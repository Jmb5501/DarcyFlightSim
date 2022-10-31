function [] = plotSimData(data)
%plotSimData plots data output byt the runSim function

f1 = figure();

subplot(3, 3, 1)
plot(data.t, convforce(data.thrust, 'N', 'lbf'));
xlabel("Time (s)");
ylabel("Thrust (lbf)");
title("Thrust vs. Time");

subplot(3, 3, 2)
plot(data.t, convpres(data.P_c*1000, 'pa', 'psi'));
xlabel("Time (s)");
ylabel("Pressure (psi)");
title("Chamber Pressure vs. Time");

subplot(3, 3, 3)
plot(data.t, data.m_dot);
xlabel("Time (s)");
ylabel("m-dot (kg/s)");
title("m-dot vs. Time");

subplot(3, 3, 4)
plot(data.t, data.OF);
xlabel("Time (s)");
ylabel("Oxygen-to-fuel Ratio");
title("OF Ratio vs. Time");

subplot(3, 3, 5)
plot(data.t, data.u_e);
xlabel("Time (s)");
ylabel("Exhaust Velocity (m/s)");
title("Exhaust Velocity vs. Time");

subplot(3, 3, 6)
plot(data.t, data.m_f, data.t, data.m_N2O);
legend(["Fuel Mass", "Oxidizer Mass"]);
xlabel("Time (s)");
ylabel("Mass (kg)");
title("Propellant Masses Remaining vs. Time");

subplot(3, 3, 7)
plot(data.t, data.dP_f, data.t, data.dP_N2O);
legend(["Fuel dP", "Oxidizer dP"]);
xlabel("Time (s)");
ylabel("Pressure (Pa)");
title("Fuel- and Oxidizer-Side Pressure Drops vs. Time");

f1.WindowState = 'maximized';

end