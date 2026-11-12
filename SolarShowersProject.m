%% Solar Shower Transient Thermal Model
% Single-day simulation
% Location: Yellowstone National Park

close all;
clear;
clc;

%% Design Parameters

% Pipe geometry
L = 1.8288;                 % Pipe length [m]
Do = 0.1016;                % Pipe outer diameter [m]
t_wall = 0.003;             % Pipe wall thickness [m]
Di = Do - 2*t_wall;         % Pipe inner diameter [m]

Ac = pi*(Do^2 - Di^2)/4;    % Pipe wall cross-sectional area [m^2]
V = pi*Di^2*L/4;            % Internal water volume [m^3]

% Pipe properties (PVC)
rho_pipe = 1400;            % Pipe density [kg/m^3]
Cp_pipe = 900;              % Pipe specific heat [J/(kg*K)]
m_pipe = rho_pipe*Ac*L;     % Pipe mass [kg]

% Water properties
rho_water = 1000;           % Water density [kg/m^3]
Cp_water = 4186;            % Water specific heat [J/(kg*K)]
m_water = rho_water*V;      % Water mass [kg]

% Thermal capacitance
C_total = m_water*Cp_water + m_pipe*Cp_pipe;

% Surface areas
A_projected = L * Do;       % Approximate area exposed to sunlight
A_surface = pi * L * Do;    % Total outside surface area (Assuming sunlight hits entire surface)

%% Solar and Thermal Parameters

absorptivity = 0.9;
emissivity = 0.9;
h_conv = 10;                 % Convection coefficient [W/(m^2*K)]
sigma = 5.670374419e-8;      % Stefan-Boltzmann constant

%% Initial Conditions and Simulation Setup

Ti = 10 + 273.15;            % Initial water/pipe temperature [K]
Ta = 15 + 273.15;            % Ambient temperature [K]
T_target = 38 + 273.15;      % Target temperature [K]

% Time setup for a full year
dt = 3600;                     % 1-hour time step [s]
simulationTime = 365*24*3600;  % Full year duration [s]
time = 0:dt:simulationTime;
timeHours = time/3600;
dayOfYear = floor(timeHours / 24) + 1; % Tracks day 1 to 365

temperature = zeros(size(time));
temperature(1) = Ti;

Q_conv = zeros(size(time));
Q_rad = zeros(size(time));
Q_solar = zeros(size(time));
G = zeros(size(time));


%% Declination and Latitude Dependence 

phi = deg2rad(44.6); % Yellowstone latitude
G_sc = 1361;         % Solar constant [W/m^2]

G = zeros(size(time));
Ta_vector = zeros(size(time));

for i = 1:length(time)
    d = dayOfYear(i);
    
    % Solar declination angle [radians]
    delta = 23.45 * sin(deg2rad(360 * (284 + d) / 365)) * (pi/180);
    
    % Sunset hour angle [radians]
    cos_omega_s = -tan(phi) * tan(delta);
    cos_omega_s = max(-1, min(1, cos_omega_s)); % Clamp for polar day/night limits
    omega_s = acos(cos_omega_s); 
    
    % Day length in hours
    day_length = 2 * (omega_s * (180/pi)) / 15;
    sunrise_h = 12 - (day_length / 2);
    sunset_h = 12 + (day_length / 2);
    
    current_hour_of_day = mod(timeHours(i), 24);
    
    % Daily extraterrestrial radiation & sine-curve distribution
    if current_hour_of_day >= sunrise_h && current_hour_of_day <= sunset_h
        theta_z = acos(sin(phi)*sin(delta) + cos(phi)*cos(delta)*cos((current_hour_of_day - 12)*15*(pi/180)));
        G_extraterrestrial = G_sc * max(0, cos(theta_z));
        
        % Approximate atmospheric attenuation (clear sky model factor ~0.7)
        G(i) = 0.7 * G_extraterrestrial; 
    else
        G(i) = 0;
    end
end
% Annual average ambient temp ~ 2.5 C, amplitude ~ 15 C
T_a_mean = 2.5 + 273.15; 
T_a_amplitude = 15;

for i = 1:length(time)
    d = dayOfYear(i);
    current_hour_of_day = mod(timeHours(i), 24);
    
    % Seasonal variation + daily thermal swing
    seasonal_temp = T_a_mean - T_a_amplitude * cos(2 * pi * d / 365);
    daily_swing = 5 * sin(2 * pi * (current_hour_of_day - 9) / 24); % peak at 3 PM
    
    Ta_vector(i) = seasonal_temp + daily_swing;
end
%% Thermal Resistance Network

% PVC thermal conductivity
k_pipe = 0.19;             % W/(m*K), approximate PVC value - Confirm this

% Internal convection coefficient
h_inside = 35;            % W/(m^2*K), initial assumption

% Pipe radii
ro = Do/2;
ri = Di/2;

% Pipe-to-water thermal resistance
R_cond = log(ro/ri)/(2*pi*k_pipe*L);

A_inside = 2*pi*ri*L;

R_inside = 1/(h_inside*A_inside);

R_pipe_water = R_cond + R_inside;

%% Two-Node Transient Simulation

% Initialize pipe and water temperatures
T_pipe = zeros(size(time));
T_water = zeros(size(time));

T_pipe(1) = Ti;
T_water(1) = Ti;

T_sky = Ta - 6;

for i = 1:length(time)-1

    % Solar heat absorbed by pipe
    Q_solar(i) = absorptivity*A_projected*G(i);

    % Convection from pipe to ambient
    Q_conv(i) = h_conv*A_surface*(T_pipe(i)-Ta);

    % Radiation from pipe to sky
    h_rad = emissivity*sigma* ...
        (T_pipe(i)+T_sky)*(T_pipe(i)^2+T_sky^2);

    Q_rad(i) = h_rad*A_surface*(T_pipe(i)-T_sky);

    % Heat transfer from pipe to water
    Q_pipe_water = (T_pipe(i)-T_water(i))/R_pipe_water;

    % Pipe energy balance
    Q_net_pipe = Q_solar(i)-Q_conv(i)-Q_rad(i)-Q_pipe_water;

    % Water energy balance
    Q_net_water = Q_pipe_water;

    % Temperature updates
    T_pipe(i+1) = T_pipe(i) + ...
        Q_net_pipe*dt/(m_pipe*Cp_pipe);

    T_water(i+1) = T_water(i) + ...
        Q_net_water*dt/(m_water*Cp_water);

end

%% Weather Conditions

G_max = 800;                 % Peak solar irradiance [W/m^2]

sunrise = 6;                 % Sunrise [hours]
sunset = 18;                 % Sunset [hours]

for i = 1:length(time)

    if timeHours(i) >= sunrise && timeHours(i) <= sunset

        G(i) = G_max*sin(pi*(timeHours(i)-sunrise) ...
            /(sunset-sunrise));

    else

        G(i) = 0;

    end

end

% Annual average ambient temp ~ 2.5 C, amplitude ~ 15 C
T_a_mean = 2.5 + 273.15; 
T_a_amplitude = 15;

for i = 1:length(time)
    d = dayOfYear(i);
    current_hour_of_day = mod(timeHours(i), 24);
    
    % Seasonal variation + daily thermal swing
    seasonal_temp = T_a_mean - T_a_amplitude * cos(2 * pi * d / 365);
    daily_swing = 5 * sin(2 * pi * (current_hour_of_day - 9) / 24); % peak at 3 PM
    
    Ta_vector(i) = seasonal_temp + daily_swing;
end

T_water_freezing = 273.15;
latent_heat_fusion = 334000; % [J/kg]

% Inside simulation loop:
if T_water(i) <= T_water_freezing && m_water > 0
    % If water reaches freezing point, lock temperature at 0C 
    % until latent heat is fully extracted or added
    T_water(i+1) = T_water_freezing;
    
end

%% Results

temperature_C = T_water-273.15;

[maxTemp,maxIndex] = max(temperature_C);

fprintf('\n--- SOLAR SHOWER RESULTS ---\n');

fprintf('Maximum temperature: %.2f C\n',maxTemp);

fprintf('Time of maximum temperature: %.2f hours\n', ...
    timeHours(maxIndex));

% Determine when target temperature is reached
targetIndex = find(temperature_C >= 38,1,'first');

if ~isempty(targetIndex)

    fprintf('Target temperature reached at %.2f hours\n', ...
        timeHours(targetIndex));

    % Calculate duration above target
    aboveTarget = temperature_C >= 38;

    totalTimeAbove = sum(aboveTarget)*dt/3600;

    fprintf('Total time above target: %.2f hours\n', ...
        totalTimeAbove);

else

    fprintf('Target temperature was NOT reached.\n');

end

% Calculate time above target after sunset
afterSunset = timeHours >= sunset;

hotAfterSunset = afterSunset & (temperature_C >= 38);

timeAfterSunset = sum(hotAfterSunset)*dt/3600;

fprintf('Time above target after sunset: %.2f hours\n', ...
    timeAfterSunset);

%% Temperature Plot

figure;

plot(timeHours,temperature_C,'LineWidth',2);
hold on;

yline(38,'r--','Target Temperature');
xline(sunrise,'k:','Sunrise');
xline(sunset,'k:','Sunset');

xlabel('Time [hours]');
ylabel('Temperature [C]');
title('Solar Shower Temperature Throughout the Day');

grid on;
xlim([0 24]);

%% Solar Irradiance Plot

figure;

plot(timeHours,G,'LineWidth',2);

xlabel('Time [hours]');
ylabel('Solar Irradiance [W/m^2]');
title('Solar Irradiance Throughout the Day');

grid on;
xlim([0 24]);
