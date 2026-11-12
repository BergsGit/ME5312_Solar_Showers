%% SOLAR SHOWER TRANSIENT THERMAL MODEL (WITH NREL DATA)
% Yellowstone National Park
% Full-year simulation using NREL TMY CSV Data

close all;
clear;
clc;

%% 1. PIPE GEOMETRY

L = 1.8288;                 % Pipe length [m]
Do = 0.1016;                % Outer diameter [m]
t_wall = 0.003;             % Pipe wall thickness [m]

Di = Do - 2*t_wall;         % Inner diameter [m]
ro = Do/2;                  % Outer radius [m]
ri = Di/2;                  % Inner radius [m]

% PVC pipe cross-sectional area
A_pipe_cross = pi*(ro^2-ri^2);

% Internal water volume
V_water = pi*Di^2*L/4;


%% 2. MATERIAL PROPERTIES

% PVC
rho_pipe = 1400;            % Density [kg/m^3]
Cp_pipe = 900;              % Specific heat [J/(kg*K)]
k_pipe = 0.19;              % Thermal conductivity [W/(m*K)]

% Water
rho_water = 1000;           % Density [kg/m^3]
Cp_water = 4186;            % Specific heat [J/(kg*K)]

% Masses & Thermal capacitances
m_pipe = rho_pipe*A_pipe_cross*L;
m_water = rho_water*V_water;
C_pipe = m_pipe*Cp_pipe;
C_water = m_water*Cp_water;


%% 3. SURFACE AREAS

A_outer = pi*Do*L;          % Outer cylindrical pipe surface [m^2]
A_inner = pi*Di*L;          % Inner cylindrical pipe surface [m^2]
A_projected = L*Do;         % Projected area for solar absorption [m^2]
A_cap_one = pi*Di^2/4;      % Area of one circular end cap [m^2]
A_caps = 2*A_cap_one;       % Total area of both end caps [m^2]


%% 4. THERMAL PROPERTIES

absorptivity = 0.90;        % Solar absorptivity
emissivity = 0.90;          % Thermal emissivity
h_inside = 35;              % Internal convection coefficient [W/(m^2*K)]
sigma = 5.670374419e-8;     % Stefan-Boltzmann constant [W/(m^2*K^4)]


%% 5. INITIAL & TARGET CONDITIONS

Ti = 10 + 273.15;           % Initial water and pipe temperature [K]
T_target = 38 + 273.15;     % Target water temperature [K]


%% 6. SIMULATION SETUP

dt = 60;                    % Reduced time step to 1 minute [s] for numerical stability
simulationTime = 365*24*3600; % Total simulation duration [s]
time = 0:dt:simulationTime;   % Time vector [s]
timeHours = time/3600;      % Time vector [hours]
N = length(time);           % Number of steps


%% 7. LOAD & INTERPOLATE NREL TMY WEATHER DATA

filename = 'Yellowstone_Data_tmy.csv';

% Explicitly specify header row at line 3 and data starting at line 4
opts = spreadsheetImportOptions('NumVariables', 10); % or generic delimited options
opts = detectImportOptions(filename, 'FileType', 'text');
opts.VariableNamingRule = 'preserve';

% Force MATLAB to read row 3 as the header names and row 4 onwards as data
nrel_data = readtable(filename, 'FileType', 'text', 'Delimiter', ',', 'HeaderLines', 2);

% Extract hourly NREL columns safely by name or column index
nrel_hours = (0:8759)' * 3600;              % 8760 hours in seconds
nrel_temp  = nrel_data.Temperature + 273.15; % Convert C to K
nrel_ghi   = nrel_data.GHI;                  % Global Horizontal Irradiance [W/m^2]

% Handle 'Wind Speed' column name safely (MATLAB converts spaces to underscores or keeps them)
if any(strcmp(nrel_data.Properties.VariableNames, 'Wind Speed'))
    nrel_wind = nrel_data.("Wind Speed");
elseif any(strcmp(nrel_data.Properties.VariableNames, 'WindSpeed'))
    nrel_wind = nrel_data.WindSpeed;
else
    nrel_wind = nrel_data{:, 10}; % Fallback to 10th column
end

% Interpolate hourly NREL data to match 5-min simulation time steps
Ta_vector   = interp1(nrel_hours, nrel_temp, time(:), 'linear', 'extrap')';
G_ghi       = interp1(nrel_hours, nrel_ghi,  time(:), 'linear', 'extrap')';
wind_vector = interp1(nrel_hours, nrel_wind, time(:), 'linear', 'extrap')';

% Ensure no negative solar values at night
G_ghi(G_ghi < 0) = 0;
G = G_ghi;

% Approximate sky temperature
Tsky_vector = Ta_vector - 6;
Tsky_vector(Tsky_vector < 1) = 1;


%% 8. THERMAL RESISTANCES

% PVC cylindrical wall conduction resistance
R_pipe_cond = log(ro/ri)/(2*pi*k_pipe*L);
R_pipe_inside = 1/(h_inside*A_inner);
R_pipe_water = R_pipe_cond + R_pipe_inside;

% End cap resistances
t_cap = t_wall;
R_cap_cond = t_cap/(k_pipe*A_cap_one);
R_caps_cond = R_cap_cond/2;


%% 9. INITIALIZE ARRAYS

T_pipe = zeros(1, N);
T_water = zeros(1, N);
Q_solar = zeros(1, N);
Q_conv = zeros(1, N);
Q_rad = zeros(1, N);
Q_pipe_water = zeros(1, N);
Q_cap = zeros(1, N);

T_pipe(1) = Ti;
T_water(1) = Ti;


%% 10. TRANSIENT THERMAL SIMULATION LOOP

for i = 1:N-1
    Tp = T_pipe(i);
    Tw = T_water(i);
    Ta = Ta_vector(i);
    Tsky = Tsky_vector(i);
    
    % Dynamic external convection coefficient based on wind speed
    V_wind = wind_vector(i);
    h_conv = max(5.0, 5.7 + 3.8 * V_wind);

    % Solar Heat Absorption
    Q_solar(i) = absorptivity * A_projected * G(i);

    % Outside Pipe Convection
    Q_conv(i) = h_conv * A_outer * (Tp - Ta);

    % Outside Pipe Radiation
    h_rad_pipe = emissivity * sigma * (Tp + Tsky) * (Tp^2 + Tsky^2);
    Q_rad(i) = h_rad_pipe * A_outer * (Tp - Tsky);

    % Pipe to Water Heat Transfer
    Q_pipe_water(i) = (Tp - Tw) / R_pipe_water;

    % End Cap Heat Loss
    T_cap = (Tw + Ta) / 2;
    h_rad_cap = emissivity * sigma * (T_cap + Tsky) * (T_cap^2 + Tsky^2);
    h_cap_external = h_conv + h_rad_cap;
    R_cap_external = 1 / (h_cap_external * A_caps);
    R_cap_total = R_caps_cond + R_cap_external;
    Q_cap(i) = (Tw - Ta) / R_cap_total;

    % Net Energy Balances
    Q_net_pipe = Q_solar(i) - Q_conv(i) - Q_rad(i) - Q_pipe_water(i);
    Q_net_water = Q_pipe_water(i) - Q_cap(i);

    % Temperature Updates
    T_pipe(i+1) = Tp + Q_net_pipe * dt / C_pipe;
    T_water(i+1) = Tw + Q_net_water * dt / C_water;
end


%% 11. CONVERT TEMPERATURES TO CELSIUS

T_pipe_C = T_pipe - 273.15;
T_water_C = T_water - 273.15;
Ta_C = Ta_vector - 273.15;

% Calculate day of year vector for monthly plots
dayOfYear = min(floor(timeHours/24) + 1, 365);


%% 12. RESULTS & SUMMARY OUTPUT

[maxTemp, maxIndex] = max(T_water_C);
fprintf('\n====================================\n');
fprintf('   SOLAR SHOWER RESULTS (NREL)\n');
fprintf('====================================\n');
fprintf('Maximum water temperature: %.2f C\n', maxTemp);
fprintf('Time of max water temp: Hour %.1f (Day %.1f)\n', timeHours(maxIndex), timeHours(maxIndex)/24);


%% 13. MONTHLY WATER TEMPERATURE SUBPLOTS (12-MONTH GRID)

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
               'July', 'August', 'September', 'October', 'November', 'December'};

start_days = [1, cumsum(days_in_month(1:end-1)) + 1];
end_days = cumsum(days_in_month);

figure('Name', 'Monthly Water Temperature Profiles', 'NumberTitle', 'off');

for m = 1:12
    % Mask for the specific month
    month_mask = (dayOfYear >= start_days(m)) & (dayOfYear <= end_days(m));
    month_days = (timeHours(month_mask) / 24) - start_days(m) + 1;
    
    subplot(4, 3, m);
    plot(month_days, T_water_C(month_mask), 'b-', 'LineWidth', 1.2);
    hold on;
    yline(37, 'r--', 'LineWidth', 1); % Target temp line
    
    title(month_names{m}, 'FontSize', 10);
    xlabel('Day of Month', 'FontSize', 8);
    ylabel('Temp [°C]', 'FontSize', 8);
    grid on;
    xlim([1, days_in_month(m)]);
    ylim([-40, 60]);
end

sgtitle('Water Temperature Profiles by Month (Yellowstone NREL Data)', 'FontSize', 12, 'FontWeight', 'bold');

%% 14. DETAILED MONTHLY METRIC ANALYSIS

days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
month_names = {'January', 'February', 'March', 'April', 'May', 'June', ...
               'July', 'August', 'September', 'October', 'November', 'December'};

start_days = [1, cumsum(days_in_month(1:end-1)) + 1];
end_days = cumsum(days_in_month);

T_target = 38; % Target temperature [°C]

fprintf('\n======================================================================================\n');
fprintf(' %-10s | %-11s | %-18s | %-18s | %-16s \n', ...
    'Month', 'Target Met?', 'Avg Heating Time', 'Max Daily Window', 'Retention Hours');
fprintf('======================================================================================\n');

for m = 1:12
    % Isolate month indices
    m_mask = (dayOfYear >= start_days(m)) & (dayOfYear <= end_days(m));
    m_tw = T_water_C(m_mask);
    m_hours = timeHours(m_mask);
    m_days = dayOfYear(m_mask);
    
    % Check if target is achieved at all during the month
    target_reached = max(m_tw) >= T_target;
    
    if target_reached
        % Track statistics on days where target was met
        heating_times = [];
        window_durations = [];
        
        for d = start_days(m):end_days(m)
            d_mask = (m_days == d);
            d_tw = m_tw(d_mask);
            d_h = m_hours(d_mask);
            
            if max(d_tw) >= T_target
                % Index when water first reaches 38°C
                idx_start = find(d_tw >= T_target, 1, 'first');
                % Index when water drops back below 38°C
                idx_end   = find(d_tw >= T_target, 1, 'last');
                
                % Heating time measured from 6:00 AM sunrise baseline
                sunrise_hour = (d - 1) * 24 + 6;
                time_to_heat = d_h(idx_start) - sunrise_hour;
                if time_to_heat > 0
                    heating_times = [heating_times, time_to_heat];
                end
                
                % Total continuous hours above 38°C on this day
                duration = (d_h(idx_end) - d_h(idx_start));
                window_durations = [window_durations, duration];
            end
        end
        
        avg_heat_time = mean(heating_times);
        max_duration = max(window_durations);
        total_month_hours = sum(m_tw >= T_target) * (dt / 3600);
        
        fprintf(' %-10s | %-12s | %-18.1f hrs | %-16.1f hrs | %-16.1f hrs \n', ...
            month_names{m}, 'YES', avg_heat_time, max_duration, total_month_hours);
    else
        fprintf(' %-10s | %-12s | %-18s | %-16s | %-16s \n', ...
            month_names{m}, 'NO', 'N/A (Too Cold)', '0.0 hrs', '0.0 hrs');
    end
end
fprintf('========================================================================================\n');