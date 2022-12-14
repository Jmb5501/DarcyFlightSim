function [params] = engineSetup()
%%%%% Configure Parameters for engineSim.m %%%%%
% Purpose: getParams Returns engine parameters for engine as "params" structure
% Original Author: Paul Sperling
% Last Modified: John Beavers
% Date Updated: 10/19/22

%% Sim Parameters
% Sim parameters are the parameters that control the simulator, like
% timestep

params.timestep = 0.1; % seconds
params.efficiency = 0.83; % Efficiency, in terms of thrust vs. ideal
params.alt = 0; % Altitude in meters, informs atmospheric pressure
params.temp = 25; % Degrees Celsius, temperature of nitrous (approx ambient temp)

%% Engine Parameters

params.fuel = "IPA"; % Stores the fuel that we are using; options are "IPA", "Ethanol,"...
[params.CDA_f, params.CDA_fI] = getFCDAs(); % Fuel CDA of feedsystem, Fuel injector CDA, calculated by other function
[params.CDA_o, params.CDA_oI] = getOCDAs(); % Oxidizer CDA of feedsystem, oxidizer injector CDA, calculated by other function
throatDiam = 0.875; % Inches, throat diameter
params.A_throat = pi()*(convlength(throatDiam, 'in', 'm')/2)^2; % Throat area in m
params.Rho_f = 786; % kg/m^3, fuel density
params.fTankPDiff = 103; % kPa, difference in tank pressure between oxidizer and fuel due to piston friction


%% Tank Params
%Calculates tank volumes in terms of lengths and diameters

concentricTanks = true; % Defines whether the fuel tank is inside the oxidizer tank

if(concentricTanks)
    oTankID = 3.75; % Inches, oxdizer tank ID
    oTankL = 72; % Inches, oxidizer tank length
    oTankUllL = 3; % Inches, length of unfilled region at top of oxidizer tank (CHECK D2 CAD)
    
    fTankOD = 1.75; % Inches, fuel tank OD
    fTankID = 1.68; % Inches, fuel tank ID
    fTankL = 70; % Inches, fuel tank length
    fPistonV = 2.0534; % Inches cubed, fuel piston displacement volume;
    fTankUllL = 3; % Inches, Length of region in top of fuel tank not filled with fuel (Fuel "ullage" length)
    
    fV = (((fTankID/2)^2)*pi()*(fTankL-fTankUllL))-fPistonV; % Inches Cubed, Fuel volume
    fTankDV = fV + pi()*(((fTankOD/2)^2)-((fTankID/2)^2))*fTankL + fPistonV; % Inches Cubed, Fuel Tank Displacement Volume
    oV = (oTankL-oTankUllL)*pi()*((oTankID/2)^2) - fTankDV; % Inches Cubes, Oxidizer Volume
    
    params.V_f = convlength(fV^(1/3), 'in', 'm')^3; % meters cubed, fuel volume
    params.V_o = convlength(oV^(1/3), 'in', 'm')^3; % meters cubed, oxidizer volume
else
    oTankID = 2.875; % Inches, oxdizer tank ID
    oTankL = 72; % Inches, oxidizer tank length
    oTankUllL = 3; % Inches, length of unfilled region at top of oxidizer tank

    fTankID = 1.375; % Inches, fuel tank ID
    fTankL = 70; % Inches, fuel tank length
    fPistonV = 2; % Inches cubed, fuel piston displacement volume;
    fTankUllL = 3; % Inches, Length of region in top of fuel tank not filled with fuel (Fuel "ullage" length)

    fV = (((fTankID/2)^2)*pi()*(fTankL-fTankUllL))-fPistonV; % Inches Cubed, Fuel volume
    oV = (oTankL - oTankIllL)*pi()*((oTankID/2)^2); % Inches Cubes, Oxidizer Volume

    params.V_f = convlength(fV^(1/3), 'in', 'm')^3; % meters cubed, fuel volume
    params.V_o = convlength(oV^(1/3), 'in', 'm')^3; % meters cubed, oxidizer volume
end

end

%%% HELPER FUNCTIONS
function [CDA_feed_f,CDA_inj_f] = getFCDAs()
%getFCDAs calculates and returns the fuel side CDA equivalent for the
%feedsystem, and the fuel side injector equivalent CDA
CDA_feed_f = 0.3*0.00000672; % Feedsystem CDA for nitrous side
CDA_inj_f = 5;
end

function [CDA_feed_N2O,CDA_inj_N20] = getOCDAs()
%getOCDAs calculates and returns the oxidizer side CDA equivalent for the
%feedsystem, and the oxidizer side injector equivalent CDA
CDA_feed_N2O = 0.25*0.00002462; % Feedsystem CDA for nitrous side
CDA_inj_N20 = 5;
end