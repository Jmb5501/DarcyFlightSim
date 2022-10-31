function [data, avgThrust, maxThrust, impulse, isp] = runEngineSim(params, verbose)
%%%%% Engine Model Darcy II%%%%%
% Purpose: runSim runs the simulator on a given set of engine parameters
% Original Author: Paul Sperling
% Last Modified: John Beavers
% Date Updated: 10/19/22


% Verbose (bool) tells the function whether or not it should display parameters
% and write verbose data to the console.
% params = getParams;
% verbose = true;

%% Setting initial parameters
t(1) = 0; % seconds, time
[P_N2O(1), rho_N2O] = getN2OPrho(params.temp); % kPa, kg/m^3, initial nitrous pressure and nitrous liquid density
m_N2O(1) = params.V_o*rho_N2O; % kg, mass of nitrous in tank initially
m_f(1) = params.V_f*params.Rho_f; % kg, mass of fuel in tank initially
P_f(1) = P_N2O(1)-params.fTankPDiff; % Initial fuel tank pressure
R_s(1) = 400; % j/kg*K, specific gas constant HARDCODED FOR NOW BUT SHOULD PROBABLY NOT BE HARDCODED
P_c(1:2) = [convpres(55, 'psi', 'Pa')/1000, convpres(500, 'psi', 'Pa')/1000]; % kPa, chamber pressure
gam_0(1:2) = [1.26, 1.2486]; % ul, ratio of specific heats for chamber must be hard-coded to start burn
gam_Th(1:2) = [1.2, 1.2]; % ul, ratio of specific heats for throat ""
gam_e(1:2) = [1.13, 1.272]; % ul, ratio of specific heats for exit ""
dP_N2O(1) = P_N2O(1) - P_c(1); % kPa, dP across nitrous side of feedsystem
dP_f(1) = P_f(1) - P_c(1); % kPa, dP across fuel side of feedsystem
m_dot_N2O(1) = params.CDA_o * sqrt(2*rho_N2O*dP_N2O(1)*1000); % kg/s, m-dot of oxidizer
m_dot_f(1) = params.CDA_f * sqrt(2*params.Rho_f*dP_f(1)*1000); % kg/s, m-dot of fuel
m_dot(1) = m_dot_f(1) + m_dot_N2O(1); % kg/s, combined m-dot of fuel and oxidizer
m_step_N2O(1) = m_dot_N2O(1)*params.timestep; % kg, change in oxidizer mass over iteration
m_step_f(1) = m_dot_f(1) * params.timestep; % kg, change in oxidizer mass over iteration
OF(1) = m_dot_N2O/m_dot_f; % Oxygen-to-fuel mass ratio
T_0(1) = 273.12 + params.temp; % Kelvins, chamber temperature
u_e(1) = 400; % m/s, exhaust velocity
thrust(1) = m_dot(1)*u_e(1); % Newtons, thrust
I_step(1) = thrust(1)*params.timestep; % newton-seconds, impulse delivered per iteration
I_tot(1) = 0; % newton-seconds, impulse delivered so far
liquidN2O(1) = 1; % boolean, there is still liquid N2O

%% Running iterative loop for engine burn
i = 2;
while ((m_N2O(i-1) > 0) && (m_f(i-1) > 0)) && (t(i-1) < 100)
    t(i) = t(i-1) + params.timestep; % seconds, update time based on timestep
    m_N2O(i) = m_N2O(i-1) - m_step_N2O(i-1); % kg, update nitrous tank mass based on amount used last iteration
    m_f(i) = m_f(i-1) - m_step_f(i-1); % kg, update fuel tank mass based on amount used last iteration
    if liquidN2O(i-1) == 1
        P_N2O(i) = P_N2O(1)*(0.7 + (m_N2O(i-1)/m_N2O(1))*0.3); % kPa, update nitrous tank pressure based on amount of nitrous left (liquid)
    else 
        P_N2O(i) = P_N2O(i-1)*exp(-0.5*(params.timestep)); % kPa, update nitrous tank pressure when it's gaseous
    end
    P_f(i) = P_N2O(i) - params.fTankPDiff; % kPa, update fuel tank pressure based on nitrous tank pressure
    if i > 2 % THE FIRST TWO ITERATIONS OF THESE ARE HARDCODED
        gam_e(i) = lookupGam_e(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel); % ratio of specific heats, exit, based on lookup table
        gam_Th(i) = lookupGam_Th(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel); % ratio of specific heats, throat, based on lookup table
        gam_0(i) = lookupGam_0(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel); % ratio of specific heats, chamber, based on lookup table
    end
    if i > 1 % THE FIRST ITERATION OF THESE ARE HARDCODED
        T_0(i) = lookupT_0(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel); % Kelvins, chamber temperature, based on lookup table
        u_e(i) = lookupv_e(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel); % m/s, exit velocity, based on lookup table
    end
    R_s(i) = lookupR(convpres(P_c(i-1)*1000, 'Pa', 'psi'), OF(i-1), params.fuel)*1000; % j/kg*K, specific gas constant, ONLY CODED PROPERLY FOR IPA AT THE MOMENT
%     R_s(i) = 400; % hardcoded for a sec
    if i>2
        P_c(i) = m_dot(i-1)*sqrt(T_0(i-1))/(params.A_throat*sqrt((gam_0(i-1)/R_s(i-1))*((2/(gam_0(i-1)+1))^((gam_0(i-1)+1)/(gam_0(i-1)-1)))))/1000; % kPa, chamber pressure derived from unknown equation that apparently is a thing
    end
    dP_N2O(i) = P_N2O(i) - P_c(i); % kPa, dP across nitrous side of feedsystem
    dP_f(i) = P_f(i) - P_c(i); % kPa, dP across fuel side of feedsystem
    m_dot_N2O(i) = params.CDA_o * sqrt(2*rho_N2O*dP_N2O(i)*1000); % kg/s, m-dot of oxidizer
    m_dot_f(i) = params.CDA_f * sqrt(2*params.Rho_f*dP_f(i)*1000); % kg/s, m-dot of fuel FIX IMAGINARY NUMBER ISSUE ON THIS LINE
    m_dot(i) = m_dot_f(i) + m_dot_N2O(i); % kg/s, combined m-dot of fuel and oxidizer
    m_step_N2O(i) = m_dot_N2O(i)*params.timestep; % kg, change in oxidizer mass over iteration
    m_step_f(i) = m_dot_f(i) * params.timestep; % kg, change in oxidizer mass over iteration
    OF(i) = m_dot_N2O(i)/m_dot_f(i); % Oxygen-to-fuel mass ratio
    thrust(i) = m_dot(i)*u_e(i);
    I_step(i) = thrust(i)*params.timestep; % newton-seconds, impulse delivered per iteration
    I_tot(i) = I_tot(i-1) + I_step(i);
    if m_N2O(i) > 0.1*m_N2O(1)
        liquidN2O(i) = 1;
    else
        liquidN2O(i) = 0;
    end
    i = i+1; % iterate
end

%% Storing data and outputting

data.t = t; 
data.m_N2O = m_N2O; 
data.m_f = m_f;
data.P_N2O = P_N2O; 
data.P_f = P_f; 
data.gam_e = gam_e; 
data.gam_Th = gam_Th;
data.gam_0 = gam_0;
data.T_0 = T_0;
data.u_e = u_e; 
data.R_s = R_s;
data.P_c = P_c;
data.dP_N2O = dP_N2O;
data.dP_f = dP_f;
data.m_dot_N2O = m_dot_N2O;
data.m_dot_f = m_dot_f;
data.m_dot = m_dot;
data.m_step_N2O = m_step_N2O; 
data.m_step_f = m_step_f;
data.OF = OF;
data.thrust = thrust;
data.I_step = I_step;
data.I_tot = I_tot; 
data.liquidN2O = liquidN2O;

avgThrust = mean(thrust); % N, Calculate and return average thrust
maxThrust = max(thrust); % N, Calculate and return max thrust
impulse = I_tot(length(I_tot)); % N-s, calculate and return total impulse;
isp = impulse/((m_N2O(1) + m_f(1))*9.8065); % Specific impulse, seconds

%% Plot data if verbose

if(verbose)
   plotSimData(data);
end

end


%% HELPER FUNCTIONS
function [pressure,rho] = getN2OPrho(T)
%getN2OP returns the nitrous pressure and density based on nitrous
%temperature (approximated as ambient)
% MUST INTAKE TEMPERATURE IN CELSIUS

load nitrousData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
while(lookupN2OTemp(i) < T)
    i = i+1;
end
i = i-1;
diff = T - lookupN2OTemp(i); 
indexT = i + (diff/(lookupN2OTemp(i+1) - lookupN2OTemp(i)));

pPpi = lookupN2OP(i+1) - lookupN2OP(i); % Partial of pressure with respect to i
pressure = lookupN2OP(i) + pPpi*(indexT - i); % kPa, nitrous pressure

prhopi = lookupN2Orho(i+1) - lookupN2Orho(i); % Partial of rho with respect to i
rho = lookupN2Orho(i) + prhopi*(indexT - i); % kg/m^3, density of liquid nitrous (for initial mass calculations)
end

function [gam_0] = lookupGam_0(P_c, OF, fuel)
%lookupGam_0 looks up the initial (chamber) ratio of specific heats of the
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    gam_0 = IPAgam_0(i, j);
    pGam_pi = IPAgam_0(i+1, j) - IPAgam_0(i, j);
    pGam_pj = IPAgam_0(i, j+1) - IPAgam_0(i, j);
    gam_0 = gam_0 + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    gam_0 = Ethanolgam_0(i, j);
    pGam_pi = Ethanolgam_0(i+1, j) - Ethanolgam_0(i, j);
    pGam_pj = Ethanolgam_0(i, j+1) - Ethanolgam_0(i, j);
    gam_0 = gam_0 + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Methanol'
    gam_0 = Methanolgam_0(i, j);
    pGam_pi = Methanolgam_0(i+1, j) - Methanolgam_0(i, j);
    pGam_pj = Methanolgam_0(i, j+1) - Methanolgam_0(i, j);
    gam_0 = gam_0 + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

function [gam_e] = lookupGam_e(P_c, OF, fuel)
%lookupGam_e looks up the exhaust ratio of specific heats of the
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    gam_e = IPAgam_e(i, j);
    pGam_pi = IPAgam_e(i+1, j) - IPAgam_e(i, j);
    pGam_pj = IPAgam_e(i, j+1) - IPAgam_e(i, j);
    gam_e = gam_e + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    gam_e = Ethanolgam_e(i, j);
    pGam_pi = Ethanolgam_e(i+1, j) - Ethanolgam_e(i, j);
    pGam_pj = Ethanolgam_e(i, j+1) - Ethanolgam_e(i, j);
    gam_e = gam_e + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Methanol'
    gam_e = Methanolgam_e(i, j);
    pGam_pi = Methanolgam_e(i+1, j) - Methanolgam_e(i, j);
    pGam_pj = Methanolgam_e(i, j+1) - Methanolgam_e(i, j);
    gam_e = gam_e + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

function [gam_Th] = lookupGam_Th(P_c, OF, fuel)
%lookupgam_Th looks up the throat ratio of specific heats of the
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    gam_Th = IPAgam_Th(i, j);
    pGam_pi = IPAgam_Th(i+1, j) - IPAgam_Th(i, j);
    pGam_pj = IPAgam_Th(i, j+1) - IPAgam_Th(i, j);
    gam_Th = gam_Th + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    gam_Th = Ethanolgam_Th(i, j);
    pGam_pi = Ethanolgam_Th(i+1, j) - Ethanolgam_Th(i, j);
    pGam_pj = Ethanolgam_Th(i, j+1) - Ethanolgam_Th(i, j);
    gam_Th = gam_Th + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

elseif fuel == 'Methanol'
    gam_Th = Methanolgam_Th(i, j);
    pGam_pi = Methanolgam_Th(i+1, j) - Methanolgam_Th(i, j);
    pGam_pj = Methanolgam_Th(i, j+1) - Methanolgam_Th(i, j);
    gam_Th = gam_Th + pGam_pi*(indexP_c - i) + pGam_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

function [R] = lookupR(P_c, OF, fuel)
%lookupT_0 looks up the specific gas constant of the 
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant
%   THE VALUE USED HERE IS AT THE THROAT, NOT THE CHAMBER

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    R = R_IPA(i, j);
    pT_pi = R_IPA(i+1, j) - R_IPA(i, j);
    pT_pj = R_IPA(i, j+1) - R_IPA(i, j);
    R = R + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    R = EthanolT_0(i, j);
    pT_pi = EthanolT_0(i+1, j) - EthanolT_0(i, j);
    pT_pj = EthanolT_0(i, j+1) - EthanolT_0(i, j);
    R = R + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

elseif fuel == 'Methanol'
    R = MethanolT_0(i, j);
    pT_pi = MethanolT_0(i+1, j) - MethanolT_0(i, j);
    pT_pj = MethanolT_0(i, j+1) - MethanolT_0(i, j);
    R = R + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

function [T_0] = lookupT_0(P_c, OF, fuel)
%lookupT_0 looks up the initial (chamber) temperature of the 
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    T_0 = IPAT_0(i, j);
    pT_pi = IPAT_0(i+1, j) - IPAT_0(i, j);
    pT_pj = IPAT_0(i, j+1) - IPAT_0(i, j);
    T_0 = T_0 + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    T_0 = EthanolT_0(i, j);
    pT_pi = EthanolT_0(i+1, j) - EthanolT_0(i, j);
    pT_pj = EthanolT_0(i, j+1) - EthanolT_0(i, j);
    T_0 = T_0 + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

elseif fuel == 'Methanol'
    T_0 = MethanolT_0(i, j);
    pT_pi = MethanolT_0(i+1, j) - MethanolT_0(i, j);
    pT_pj = MethanolT_0(i, j+1) - MethanolT_0(i, j);
    T_0 = T_0 + pT_pi*(indexP_c - i) + pT_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

function [v_e] = lookupv_e(P_c, OF, fuel)
%lookupv_e looks up the initial (chamber) temperature of the 
%flow based on the chamber pressure, the O/F ratio, and the fuel
%   First, the function looks up the index of a certain pressure and
%   mixture ratio along the "X and Y" axis of an empirical table, as a
%   fraction of an index (i.e. the value can be between two ints). Then, it
%   uses these indices to look up/interpolate a value on a 3-D surface, and
%   returns the "Z" value. The data tables are different for each
%   propellant

load propData.mat;

% Find the "index" of the chamber pressure (non-int)
i=1;
if (P_c > lookupP_c(1)) && (P_c < lookupP_c(length(lookupP_c)))
    while(lookupP_c(i) < P_c)
        i = i+1;
    end
    i = i-1;
    diff = P_c - lookupP_c(i); 
    indexP_c = i + (diff/(lookupP_c(i+1) - lookupP_c(i)));
elseif (P_c < lookupP_c(1))
    i = 1;
    indexP_c = i;
else
    i = length(lookupP_c)-1;
    indexP_c = i;
end

% Find the "index" of the O/F Ratio (non-int)
j=1;
if(OF > lookupOF(1)) && (OF < lookupOF(length(lookupOF)))
    while(lookupOF(j) < OF)
        j = j+1;
    end
    j = j-1;
    diff = OF - lookupOF(j); 
    indexOF = j + (diff/(lookupOF(j+1) - lookupOF(j)));
elseif(OF < lookupOF(1))
    i=1;
    indexOF = i;
else
    i = length(lookupOF)-1;
    indexOF = i;
end

if fuel == 'IPA'
    v_e = IPAv_e(i, j);
    pv_pi = IPAv_e(i+1, j) - IPAv_e(i, j);
    pv_pj = IPAv_e(i, j+1) - IPAv_e(i, j);
    v_e = v_e + pv_pi*(indexP_c - i) + pv_pj*(indexOF - j);

elseif fuel == 'Ethanol'
    v_e = Ethanolv_e(i, j);
    pv_pi = Ethanolv_e(i+1, j) - Ethanolv_e(i, j);
    pv_pj = Ethanolv_e(i, j+1) - Ethanolv_e(i, j);
    v_e = v_e + pv_pi*(indexP_c - i) + pv_pj*(indexOF - j);

elseif fuel == 'Methanol'
    v_e = Methanolv_e(i, j);
    pv_pi = Methanolv_e(i+1, j) - Methanolv_e(i, j);
    pv_pj = Methanolv_e(i, j+1) - Methanolv_e(i, j);
    v_e = v_e + pv_pi*(indexP_c - i) + pv_pj*(indexOF - j);

else
    disp("Please input a valid fuel value (IPA, Ethanol, or Methanol")
end
end

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