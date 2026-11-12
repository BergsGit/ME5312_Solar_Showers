%% SOLAR SHOWER TRANSIENT THERMAL MODEL - REVISED
% Yellowstone, Wyoming | NSRDB Typical Meteorological Year data
%
% This version corrects the weather-file import and timing, uses DNI/DHI/GHI
% when calculating solar input, estimates solar position and sunrise, adds
% internal convection to the end-cap resistance, and improves monthly
% performance metrics.
%
% IMPORTANT ASSUMPTIONS:
% 1) One uniform temperature represents all PVC; one uniform temperature
%    represents all water (lumped-capacitance model).
% 2) The pipe is horizontal. Its axis azimuth is set below.
% 3) Internal water convection is represented by an assumed constant h_inside.
% 4) The sky temperature is approximated as ambient temperature minus 6 K. ****WE SHOULD CHANGE THIS****
% 5) Diffuse sky and ground-reflected radiation are treated as isotropic.
% 6) Weather records are hourly averages stamped at HH:30.
% 7) TMY records are a representative year assembled from source years;
%    month/day are used for solar position, not the source year.

close all;
clear;
clc;

%% 1. USER SETTINGS

filename = 'Yellowstone_Data_tmy.csv';

% Location taken from the CSV metadata
latitude_deg = 44.53;
longitude_deg = -110.82;
utc_offset_hours = -7;        % Local standard time, UTC-7

% Pipe axis azimuth, clockwise from north:
% 0 = north-south; 90 = east-west
pipe_axis_azimuth_deg = 90;

ground_reflectance = 0.20;    % Assumed ground albedo [-]

dt = 60;                      % Simulation timestep [s]
initial_temperature_C = 10;   % Initial pipe and water temperature [deg C]
target_temperature_C = 38;    % Desired water temperature [deg C]

%% 2. PIPE GEOMETRY

L = 1.8288;                   % Pipe length [m]
Do = 0.1016;                  % Outer diameter [m]
t_wall = 0.003;               % Pipe wall thickness [m]

Di = Do - 2*t_wall;
ro = Do/2;
ri = Di/2;

if Di <= 0
    error('Pipe wall thickness must be less than half the outer diameter.');
end

A_pipe_cross = pi*(ro^2 - ri^2);
V_water = pi*Di^2*L/4;

A_outer = pi*Do*L;            % Outer cylindrical area [m^2]
A_inner = pi*Di*L;            % Inner cylindrical area [m^2]
A_projected_GHI = Do*L;       % Approximate projected area for fallback only
A_caps = 2*pi*Di^2/4;         % Both internal end-cap areas [m^2]

%% 3. MATERIAL PROPERTIES

rho_pipe = 1400;              % PVC density [kg/m^3]
Cp_pipe = 900;                % PVC specific heat [J/(kg*K)]
k_pipe = 0.19;                % PVC thermal conductivity [W/(m*K)]

rho_water = 1000;             % Water density [kg/m^3]
Cp_water = 4186;              % Water specific heat [J/(kg*K)]

m_pipe = rho_pipe*A_pipe_cross*L;
m_water = rho_water*V_water;

C_pipe = m_pipe*Cp_pipe;      % Pipe thermal capacitance [J/K]
C_water = m_water*Cp_water;   % Water thermal capacitance [J/K]

%% 4. HEAT-TRANSFER PARAMETERS

absorptivity = 0.90;
emissivity = 0.90;
h_inside = 35;                % ASSUMED internal convection [W/(m^2*K)]
sigma = 5.670374419e-8;       % Stefan-Boltzmann constant [W/(m^2*K^4)]

% Empirical external wind convection approximation.
% CSV wind speed is in m/s. Validate this correlation for your application.
h_wind = @(V) max(5.0, 5.7 + 3.8*max(V,0));

%% 5. INITIAL CONDITIONS

Ti_K = initial_temperature_C + 273.15;
target_temperature_K = target_temperature_C + 273.15; %#ok<NASGU>

%% 6. READ THE NSRDB CSV FILE

if ~isfile(filename)
    error('Cannot find %s. Put it in the MATLAB Current Folder.', filename);
end

% The CSV has:
% line 1 = metadata field names
% line 2 = metadata values
% line 3 = actual weather-data column names
% line 4 onward = 8760 hourly weather records
opts = delimitedTextImportOptions('NumVariables',10);
opts.DataLines = [4 Inf];
opts.VariableNamesLine = 3;
opts.Delimiter = ',';
opts.VariableNamingRule = 'preserve';

weather = readtable(filename,opts);

requiredNames = {'Year','Month','Day','Hour','Minute', ...
                 'Temperature','DHI','DNI','GHI','Wind Speed'};

if ~all(ismember(requiredNames,weather.Properties.VariableNames))
    disp(weather.Properties.VariableNames');
    error('CSV headers do not match the expected NSRDB column names.');
end

if height(weather) ~= 8760
    warning('Expected 8760 records, but imported %d.',height(weather));
end

% Extract columns. NSRDB fields in this file:
% Temperature [deg C], DHI/DNI/GHI [W/m^2], Wind Speed [m/s].
monthData = double(weather.Month);
dayData = double(weather.Day);
hourData = double(weather.Hour);
minuteData = double(weather.("Minute"));

Ta_hourly_C = double(weather.Temperature);
DHI_hourly = max(double(weather.DHI),0);
DNI_hourly = max(double(weather.DNI),0);
GHI_hourly = max(double(weather.GHI),0);
wind_hourly = max(double(weather.("Wind Speed")),0);

% Replace missing numeric values, if any.
Ta_hourly_C = fillmissing(Ta_hourly_C,'linear','EndValues','nearest');
DHI_hourly = fillmissing(DHI_hourly,'constant',0);
DNI_hourly = fillmissing(DNI_hourly,'constant',0);
GHI_hourly = fillmissing(GHI_hourly,'constant',0);
wind_hourly = fillmissing(wind_hourly,'constant',0);

% Confirm expected hourly sequence and half-hour timestamp.
if any(hourData ~= (0:8759)' - 24*floor((0:8759)'/24)) || ...
        any(minuteData ~= 30)
    warning(['The weather timestamps are not the expected sequence of ', ...
             'hourly records at HH:30. Check the imported time columns.']);
end

%% 7. TIME VECTOR AND WEATHER INTERPOLATION

nWeather = height(weather);
simulationTime = nWeather*3600;     % 8760 hours [s]
time = 0:dt:simulationTime;
timeHours = time/3600;
N = numel(time);

% The CSV measurements are centered at 0.5, 1.5, ..., 8759.5 hours.
% Add periodic points either side so the TMY year wraps smoothly and no
% extrapolation is performed at the start/end of the simulation.
weatherHours = (0:nWeather-1)' + 0.5;
weatherHoursExtended = [-0.5; weatherHours; nWeather+0.5];

Ta_extended = [Ta_hourly_C(end); Ta_hourly_C; Ta_hourly_C(1)];
DHI_extended = [DHI_hourly(end); DHI_hourly; DHI_hourly(1)];
DNI_extended = [DNI_hourly(end); DNI_hourly; DNI_hourly(1)];
GHI_extended = [GHI_hourly(end); GHI_hourly; GHI_hourly(1)];
wind_extended = [wind_hourly(end); wind_hourly; wind_hourly(1)];

Ta_C = interp1(weatherHoursExtended,Ta_extended,timeHours,'linear');
DHI = interp1(weatherHoursExtended,DHI_extended,timeHours,'linear');
DNI = interp1(weatherHoursExtended,DNI_extended,timeHours,'linear');
GHI = interp1(weatherHoursExtended,GHI_extended,timeHours,'linear');
wind = interp1(weatherHoursExtended,wind_extended,timeHours,'linear');

Ta_K = Ta_C + 273.15;

%% 8. SOLAR POSITION, SUNRISE, AND COLLECTOR IRRADIANCE

% The source year varies by month in TMY data. Use a non-leap typical year
% based on month/day to calculate solar position.
% Solar position at every simulation time.
[solarElevation_deg,solarAzimuth_deg,sunriseHourByDay] = ...
    solarPositionTMY(timeHours,latitude_deg,longitude_deg,utc_offset_hours);

% Direct beam on a horizontal cylindrical pipe:
% projected silhouette area = D*L*sqrt(1-(sun direction dot pipe axis)^2)
sunElevation = deg2rad(solarElevation_deg);
sunAzimuth = deg2rad(solarAzimuth_deg);
pipeAxisAzimuth = deg2rad(pipe_axis_azimuth_deg);

sunAxisDot = cos(sunElevation).*cos(sunAzimuth-pipeAxisAzimuth);
sunAxisDot = min(max(sunAxisDot,-1),1);
A_beam_projected = Do*L.*sqrt(max(0,1-sunAxisDot.^2));

sunAboveHorizon = solarElevation_deg > 0;

Q_solar_direct = absorptivity .* DNI .* A_beam_projected .* ...
                 sunAboveHorizon;

% Approximate isotropic diffuse sky irradiance on the curved cylinder.
% For a vertical-facing surface, isotropic-sky diffuse irradiance is DHI/2.
Q_solar_diffuse = absorptivity .* A_outer .* 0.5 .* DHI;

% Approximate isotropic ground-reflected irradiance on the cylinder.
Q_solar_ground = absorptivity .* A_outer .* 0.5 .* ...
                 ground_reflectance .* GHI;

Q_solar_input = Q_solar_direct + Q_solar_diffuse + Q_solar_ground;

%% 9. THERMAL RESISTANCES

% Pipe wall conduction + inside convection.
R_pipe_cond = log(ro/ri)/(2*pi*k_pipe*L);
R_pipe_inside = 1/(h_inside*A_inner);
R_pipe_water = R_pipe_cond + R_pipe_inside;

% End-cap resistance includes water-to-inner-cap convection,
% PVC conduction, and external convection/radiation.
R_cap_inside = 1/(h_inside*A_caps);
R_cap_cond = t_wall/(k_pipe*A_caps);

%% 10. INITIALIZE ARRAYS

T_pipe = zeros(1,N);
T_water = zeros(1,N);

Q_solar = zeros(1,N);
Q_solar_direct_store = zeros(1,N);
Q_solar_diffuse_store = zeros(1,N);
Q_solar_ground_store = zeros(1,N);
Q_conv = zeros(1,N);
Q_rad = zeros(1,N);
Q_pipe_water = zeros(1,N);
Q_cap = zeros(1,N);

T_pipe(1) = Ti_K;
T_water(1) = Ti_K;

%% 11. TRANSIENT THERMAL SIMULATION

for i = 1:N-1

    Tp = T_pipe(i);
    Tw = T_water(i);
    Ta = Ta_K(i);

    % Simplified effective sky temperature.
    Tsky = max(Ta - 6,1);

    % External convection coefficient based on wind speed.
    h_conv = h_wind(wind(i));

    % Solar energy absorbed by the pipe [W].
    Q_solar(i) = Q_solar_input(i);
    Q_solar_direct_store(i) = Q_solar_direct(i);
    Q_solar_diffuse_store(i) = Q_solar_diffuse(i);
    Q_solar_ground_store(i) = Q_solar_ground(i);

    % External convection: positive means heat leaves the pipe.
    Q_conv(i) = h_conv*A_outer*(Tp-Ta);

    % External longwave radiation: positive means heat leaves the pipe.
    h_rad_pipe = emissivity*sigma*(Tp+Tsky)*(Tp^2+Tsky^2);
    Q_rad(i) = h_rad_pipe*A_outer*(Tp-Tsky);

    % Heat transfer from pipe wall into water.
    Q_pipe_water(i) = (Tp-Tw)/R_pipe_water;

    % End-cap loss. Linearize radiation about the water/end-cap temperature.
    h_rad_cap = emissivity*sigma*(Tw+Tsky)*(Tw^2+Tsky^2);
    h_external_cap = h_conv+h_rad_cap;
    T_external_equiv = (h_conv*Ta+h_rad_cap*Tsky)/h_external_cap;

    R_cap_external = 1/(h_external_cap*A_caps);
    R_cap_total = R_cap_inside+R_cap_cond+R_cap_external;
    Q_cap(i) = (Tw-T_external_equiv)/R_cap_total;

    % Separate lumped energy balances.
    Q_net_pipe = Q_solar(i)-Q_conv(i)-Q_rad(i)-Q_pipe_water(i);
    Q_net_water = Q_pipe_water(i)-Q_cap(i);

    % Explicit Euler updates.
    T_pipe(i+1) = Tp+(Q_net_pipe*dt/C_pipe);
    T_water(i+1) = Tw+(Q_net_water*dt/C_water);
end

% Copy final sample values for complete plotting/integration arrays.
Q_solar(end) = Q_solar_input(end);
Q_solar_direct_store(end) = Q_solar_direct(end);
Q_solar_diffuse_store(end) = Q_solar_diffuse(end);
Q_solar_ground_store(end) = Q_solar_ground(end);
Q_conv(end) = Q_conv(end-1);
Q_rad(end) = Q_rad(end-1);
Q_pipe_water(end) = Q_pipe_water(end-1);
Q_cap(end) = Q_cap(end-1);

%% 12. CONVERT TEMPERATURES AND PREPARE CALENDAR ARRAYS

T_pipe_C = T_pipe-273.15;
T_water_C = T_water-273.15;

daysInMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
monthNames = {'January','February','March','April','May','June', ...
              'July','August','September','October','November','December'};

monthStartDay = [1 cumsum(daysInMonth(1:end-1))+1];
monthEndDay = cumsum(daysInMonth);

dayOfYear = min(floor(timeHours/24)+1,365);
hourOfDay = mod(timeHours,24);

%% 13. ANNUAL SUMMARY AND ENERGY TOTALS

[maxWaterTemp_C,maxIndex] = max(T_water_C);

fprintf('\n====================================================\n');
fprintf(' SOLAR SHOWER RESULTS - YELLOWSTONE NSRDB TMY\n');
fprintf('====================================================\n');
fprintf('Weather records: %d\n',nWeather);
fprintf('Solar inputs: DNI + DHI + GHI ground reflection\n');
fprintf('Maximum water temperature: %.2f deg C\n',maxWaterTemp_C);
fprintf('Time of maximum: Day %.1f, hour %.2f\n', ...
        timeHours(maxIndex)/24+1,mod(timeHours(maxIndex),24));
fprintf('Water mass: %.2f kg\n',m_water);
fprintf('PVC mass: %.2f kg\n',m_pipe);
fprintf('Water volume: %.2f L\n',V_water*1000);
fprintf('Pipe thermal capacitance: %.1f J/K\n',C_pipe);
fprintf('Water thermal capacitance: %.1f J/K\n',C_water);

E_solar_kWh = sum(Q_solar(1:end-1))*dt/3.6e6;
E_conv_kWh = sum(Q_conv(1:end-1))*dt/3.6e6;
E_rad_kWh = sum(Q_rad(1:end-1))*dt/3.6e6;
E_pipeWater_kWh = sum(Q_pipe_water(1:end-1))*dt/3.6e6;
E_cap_kWh = sum(Q_cap(1:end-1))*dt/3.6e6;

fprintf('\nAnnual integrated energy terms (kWh):\n');
fprintf('Absorbed solar energy: %.2f\n',E_solar_kWh);
fprintf('Net external convection from pipe: %.2f\n',E_conv_kWh);
fprintf('Net radiation from pipe: %.2f\n',E_rad_kWh);
fprintf('Net pipe-to-water transfer: %.2f\n',E_pipeWater_kWh);
fprintf('Net end-cap heat loss: %.2f\n',E_cap_kWh);

%% 14. MONTHLY WATER TEMPERATURE PLOTS

figure('Name','Monthly Water Temperature Profiles','NumberTitle','off');

for m = 1:12
    monthMask = dayOfYear>=monthStartDay(m) & dayOfYear<=monthEndDay(m);
    monthDay = timeHours(monthMask)/24-monthStartDay(m)+1;

    subplot(4,3,m);
    plot(monthDay,T_water_C(monthMask),'b-','LineWidth',1.1);
    hold on;
    yline(target_temperature_C,'r--','Target','LineWidth',1);
    grid on;
    xlabel('Day of month');
    ylabel('Water temperature [deg C]');
    title(monthNames{m});
    xlim([1 daysInMonth(m)]);
end

sgtitle('Predicted Solar Shower Water Temperature by Month');

%% 15. DAILY AND MONTHLY PERFORMANCE METRICS

dailyHeatingTime_hr = nan(365,1);
dailyLongestWindow_hr = zeros(365,1);
dailyTotalAboveTarget_hr = zeros(365,1);
dailyTargetReached = false(365,1);

for d = 1:365

    dayMask = timeHours >= (d-1) & timeHours < d;
    tDay = timeHours(dayMask);
    TwDay = T_water_C(dayMask);

    if isempty(tDay)
        continue;
    end

    % Find first time at/after calculated sunrise when water is at target.
    sunriseLocal = sunriseHourByDay(d);
    afterSunrise = tDay >= ((d-1)+sunriseLocal);
    targetAfterSunrise = find(afterSunrise & (TwDay>=target_temperature_C),1,'first');

    if ~isempty(targetAfterSunrise)
        dailyTargetReached(d) = true;
        dailyHeatingTime_hr(d) = ...
            tDay(targetAfterSunrise)-((d-1)+sunriseLocal);
    end

    above = TwDay>=target_temperature_C;
    dailyTotalAboveTarget_hr(d) = sum(above)*dt/3600;

    % Identify separate continuous intervals above target.
    transitions = diff([false; above(:); false]);
    runStarts = find(transitions==1);
    runEnds = find(transitions==-1)-1;

    if ~isempty(runStarts)
        runDurations = (runEnds-runStarts+1)*dt/3600;
        dailyLongestWindow_hr(d) = max(runDurations);
    end
end

fprintf('\n====================================================================================================\n');
fprintf('%-11s | %-13s | %-18s | %-18s | %-18s | %-14s\n', ...
    'Month','Successful days','Avg heat time','Longest window','Total above target','Success rate');
fprintf('====================================================================================================\n');

for m = 1:12

    monthDays = monthStartDay(m):monthEndDay(m);
    successful = dailyTargetReached(monthDays);

    nSuccessful = sum(successful);
    successRate = 100*nSuccessful/numel(monthDays);

    if nSuccessful > 0
        avgHeating = mean(dailyHeatingTime_hr(monthDays(successful)));
    else
        avgHeating = NaN;
    end

    longestWindow = max(dailyLongestWindow_hr(monthDays));
    totalAbove = sum(dailyTotalAboveTarget_hr(monthDays));

    if isnan(avgHeating)
        avgHeatingText = 'N/A';
    else
        avgHeatingText = sprintf('%.2f hr',avgHeating);
    end

    fprintf('%-11s | %-13d | %-18s | %-18.2f | %-18.2f | %-13.1f%%\n', ...
        monthNames{m},nSuccessful,avgHeatingText,longestWindow,totalAbove,successRate);
end

fprintf('====================================================================================================\n');

%% 16. PLOT ENERGY RATES

figure('Name','Solar Shower Heat Transfer Rates','NumberTitle','off');
plot(timeHours/24,Q_solar,'LineWidth',1.0);
hold on;
plot(timeHours/24,Q_conv,'LineWidth',1.0);
plot(timeHours/24,Q_rad,'LineWidth',1.0);
plot(timeHours/24,Q_pipe_water,'LineWidth',1.0);
plot(timeHours/24,Q_cap,'LineWidth',1.0);
grid on;
xlabel('Day of typical year');
ylabel('Heat transfer rate [W]');
title('Solar Shower Heat Transfer Rates');
legend('Absorbed solar','Pipe convection','Pipe radiation', ...
       'Pipe to water','End-cap loss','Location','best');

%% 17. OPTIONAL: CHECK TIME-STEP SENSITIVITY MANUALLY
% Run the model with dt = 60, then 30, then 10 seconds.
% Compare maximum water temperature and monthly heating-time metrics.
% If results change materially as dt decreases, use the smaller timestep.
fprintf('\nTimestep used: %.1f seconds. Repeat with 30 s and 10 s to check convergence.\n',dt);

%% LOCAL FUNCTION: APPROXIMATE SOLAR POSITION FOR A TMY YEAR

function [elevationDeg,azimuthDeg,sunriseHour] = ...
    solarPositionTMY(timeHours,latitudeDeg,longitudeDeg,utcOffsetHours)

    lat = deg2rad(latitudeDeg);
    dayNumber = min(floor(timeHours/24)+1,365);
    localHour = mod(timeHours,24);

    % NOAA fractional-year approximation, using a non-leap year.
    gamma = 2*pi/365.*(dayNumber-1+(localHour-12)/24);

    equationTime = 229.18.*(0.000075 + 0.001868.*cos(gamma) ...
        - 0.032077.*sin(gamma) - 0.014615.*cos(2*gamma) ...
        - 0.040849.*sin(2*gamma));

    declination = 0.006918 - 0.399912.*cos(gamma) ...
        + 0.070257.*sin(gamma) - 0.006758.*cos(2*gamma) ...
        + 0.000907.*sin(2*gamma) - 0.002697.*cos(3*gamma) ...
        + 0.00148.*sin(3*gamma);

    standardMeridian = 15*utcOffsetHours;
    timeCorrection_min = equationTime + 4*(longitudeDeg-standardMeridian);
    solarTime = localHour + timeCorrection_min/60;
    hourAngle = deg2rad(15*(solarTime-12));

    sinElevation = sin(lat).*sin(declination) + ...
        cos(lat).*cos(declination).*cos(hourAngle);
    sinElevation = min(max(sinElevation,-1),1);
    elevation = asin(sinElevation);

    % Solar azimuth clockwise from north.
    azimuth = atan2(sin(hourAngle), ...
        cos(hourAngle).*sin(lat)-tan(declination).*cos(lat)) + pi;
    azimuth = mod(azimuth,2*pi);

    elevationDeg = rad2deg(elevation);
    azimuthDeg = rad2deg(azimuth);

    % Sunrise for each day using solar altitude = 0 degrees.
    sunriseHour = nan(365,1);

    for d = 1:365
        gammaNoon = 2*pi/365*(d-1);
        EoT = 229.18*(0.000075 + 0.001868*cos(gammaNoon) ...
            - 0.032077*sin(gammaNoon) - 0.014615*cos(2*gammaNoon) ...
            - 0.040849*sin(2*gammaNoon));

        dec = 0.006918 - 0.399912*cos(gammaNoon) ...
            + 0.070257*sin(gammaNoon) - 0.006758*cos(2*gammaNoon) ...
            + 0.000907*sin(2*gammaNoon) - 0.002697*cos(3*gammaNoon) ...
            + 0.00148*sin(3*gammaNoon);

        cosH0 = -tan(lat)*tan(dec);
        cosH0 = min(max(cosH0,-1),1);
        H0 = acos(cosH0);

        solarNoonLocal = 12-(EoT+4*(longitudeDeg-standardMeridian))/60;
        sunriseHour(d) = solarNoonLocal-rad2deg(H0)/15;
    end
end
