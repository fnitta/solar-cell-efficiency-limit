% This is the MATLAB code to calculate the fundamental performance limits of
% single-junction solar cells with a realistic analysis based on the
% Tiedje-Yablonovitch model and including defect-assisted Shockley-Read-Hall
% (SRH) recombination.
%
% See the paper here: https://www.nature.com/articles/s42005-023-01447-y

% Reset MATLAB
clear; close all;
FIGURE_COUNTER = 0;

% Define constants
planck_constant = 6.626 * 10 ^ -34; % units: J * s
electron_charge = 1.602 * 10 ^ -19; % units: C
electron_mass = 9.11 * 10 ^ -31; % units: kg
planck_eV = planck_constant / electron_charge; % units J * s * eV / J
speed_of_light = 3.0 * 10 ^ 8; % units: m / s
Boltzmann_constant = 1.381 * 10 ^ -23; % units: J / K
T = 300; % units: K

%% Section 0: Setup directories and change parameters here
% This is where you define the main folder for the modeling, which should
% contain the optical data for the material (k and n data), spectrum data,
% and this script. You should also input the electrical data for the
% material here (band gap, electron and hole effective masses, and Auger
% coefficient).

% Define main directory and data output directories
Main_Directory = "/Users/frederick.nitta/Desktop/solar-cell-efficiency-limit/";
[~, ~, ~] = mkdir("Figures (.fig)");
Figures_fig_Directory = "Figures (.fig)/";
[~, ~, ~] = mkdir("Figures (.jpg)");
Figures_jpg_Directory = "Figures (.jpg)/";
[~, ~, ~] = mkdir("Output Data");
Output_Data_Directory = "Output Data/";

% Define material
Material_Name = "WS_2";

% Define material's electrical parameters here
Material_BandGap = 1.36; % units: eV
Material_me = 0.63 * electron_mass;
Material_mh = 0.84 * electron_mass;

% Calculate effective conduction and valence band density of states and
% intrinsic carrier concentration
Material_NC = 10 ^ -6 * 2 * ((2 * pi * Material_me * Boltzmann_constant * T) / ...
    (planck_constant) ^ 2) ^ (3 / 2); % units: cm ^ -3
Material_NV = 10 ^ -6 * 2 * ((2 * pi * Material_mh * Boltzmann_constant * T) / ...
    (planck_constant) ^ 2) ^ (3 / 2); % units: cm ^ -3
Material_ni = sqrt(Material_NC * Material_NV) * exp(-(Material_BandGap * ...
    electron_charge) / (2 * Boltzmann_constant * T)); % units: cm ^ -3

% Define material's optical parameters
Material_Directory = strcat(Main_Directory, "Data/Material/", Material_Name, "/");
Material_k_File = strcat(Material_Directory, Material_Name, "-k.txt");
Material_k_Data = readmatrix(Material_k_File);
Material_n_File = strcat(Material_Directory, Material_Name, "-n.txt");
Material_n_Data = readmatrix(Material_n_File);

% Calculate alpha_2 and find n at the material's band gap energy
Material_Alpha2_Original = 4 * pi * Material_k_Data(:, 2:2) ./ Material_k_Data(:, 1:1); % units: nm ^ -1
Material_BandGap_Wavelength = (planck_eV * speed_of_light / ...
    Material_BandGap) * 10 ^ 9; % units: nm
Material_n = interp1(Material_n_Data(:, 1:1), Material_n_Data(:, 2:2), ...
    Material_BandGap_Wavelength); % no units
Material_4n2 = 4 * Material_n ^ 2; % 4 * n ^ 2, also no units
Material_n2 = Material_n ^ 2; % n ^ 2, also no units

% Define the spectrum (AM 1.5 G here)
% Downloaded from https://www.nrel.gov/grid/solar-resource/spectra-am1.5.html
Spectrum_Directory = strcat(Main_Directory, "Data/Spectrum/AM 1.5 G/");
Spectrum_File = strcat(Spectrum_Directory, "ASTMG173.csv");
Spectrum_Data = readmatrix(Spectrum_File, 'NumHeaderLines', 1);

% Make sure that the material's alpha_2 data and spectrum data are as a function
% of the same wavelength domain
Minimum_Wavelength = 280.0;
Maximum_Wavelength = 4000.0;
New_Wavelength_Range = linspace(Minimum_Wavelength, Maximum_Wavelength, ...
    20 * (Maximum_Wavelength - Minimum_Wavelength) + 1).';
New_Energy_J_Range = 10 ^ 9 * planck_constant * speed_of_light./ ...
    New_Wavelength_Range; % units: (J * s * m / s) / (nm * 10 ^ -9 m / nm) = J / photon
New_Energy_eV_range = New_Energy_J_Range / electron_charge;
Material_Alpha2 = interp1(New_Wavelength_Range, Material_Alpha2_Original, ...
    New_Wavelength_Range); % no units
Spectrum = interp1(Spectrum_Data(:, 1:1), Spectrum_Data(:, 3:3), ...
    New_Wavelength_Range); % units: W * m ^ -2 * nm ^ -1

% Define recombination parameters
Material_Alpha1 = 0; % Free Carrier Absorption coefficient; units: nm ^ -1
Material_C = 10 ^ -29.97; % Auger coefficient; units: cm ^ 6 * s ^ -1
SRH_tau_range = [Inf 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 1e-9]; % units: s

% Supports up to eight tau_SRH values; if you would like to use more,
% add more colors to this matrix
Colors = [[0.5 0 0]; [0.64 0.08 0.18]; [0.85 0.33 0.10];
    [0.93 0.69 0.13]; [0 0.5 0]; [0 0.75 0.75];
    [0.00 0.45 0.74]; [0.49 0.18 0.56]];

% Define thickness range and specific thicknesses of interest to generate
% JV curves and spectral dependencies of the luminescent emission rate
L_Range = linspace(1, 1001, 1001).'; % units: nm
L_Interest = [100 1000]; % units: nm (put in ascending order)
for i = 1:length(L_Interest)
    assert(ismember(L_Interest(i), L_Range), ...
        'Your thickness of interest is not defined in your thickness range.');
end

% Time loggers to approximate time remaining
Time_Logger_Regular = 0;
Time_Logger_Interest = 0;
%% Section 1: Calculate efficiency limits with Shockley-Queisser (SQ) Model

SQ_Voltage = linspace(0, 1.5, 1.5e3).'; % Voltage domain, units: V
Black_Body_Photon_Flux = @(E_gap) (E_gap .^ 2) ./ (exp(E_gap ./ ((Boltzmann_constant ...
    / electron_charge) * T)) - 1); % units: E_gap in eV --> eV ^ 3

% Radiative recombination
RR0 = (10 ^ 3 * 10 ^ -4 * electron_charge * 2 * pi) / (speed_of_light ^ 2 ...
    * (planck_constant / electron_charge) ^ 3) * integral(Black_Body_Photon_Flux, ...
    Material_BandGap, Inf); % units: A / m ^ 2 = 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2

% Smallest wavelength (or largest energy) that will be absorbed
Material_Lower_Bound_Wavelength = New_Wavelength_Range(1);
% Largest wavelength (or smallest energy) that will be absorbed, up to the
% band gap wavelength
Material_Upper_Bound_Wavelength = Material_BandGap_Wavelength;
% Wavelength range that will be absorbed
Material_Wavelength_Domain = linspace(Material_Lower_Bound_Wavelength, ...
    Material_Upper_Bound_Wavelength, 20 * (Material_Upper_Bound_Wavelength ...
    - Material_Lower_Bound_Wavelength) + 1).';
% Convert to wavelength to photon energy
Material_Photon_Energy = 10 ^ 9 * planck_constant * speed_of_light ./ Material_Wavelength_Domain;
% Select portion of spectrum that will be absorbed
Material_Spectrum = interp1(New_Wavelength_Range, Spectrum, Material_Wavelength_Domain); % units: J / (nm * m ^ 2)

% Calculate light-generated current by integrating the portion of the
% spectrum that is absorbed
SQ_L = 10 ^ 3 * 10 ^ -4 * electron_charge * trapz(Material_Wavelength_Domain, ...
    Material_Spectrum ./ Material_Photon_Energy); % units: (photon * C) / (nm * m ^ 2) = A / m ^ 2 =
    % 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2

% Find current density from radiative recombination and light-generated
% current, using the ideal diode equation
SQ_Current_Density = RR0 * (exp((electron_charge .* SQ_Voltage) ...
    ./ (Boltzmann_constant * T))) - SQ_L;

% Calculate SQ parameters
SQ_Jsc = abs(interp1(SQ_Voltage, SQ_Current_Density, 0));
SQ_Current_Density_new = SQ_Current_Density;
[SQ_Current_Density_new, index] = unique(SQ_Current_Density_new);
SQ_Voc = interp1(SQ_Current_Density_new, SQ_Voltage(index), 0); % V
SQ_FF = abs(min(SQ_Current_Density .* SQ_Voltage)) / abs(SQ_Voc * SQ_Jsc);
SQ_Pmax = abs(SQ_Jsc * SQ_Voc * SQ_FF); % mW / cm ^ 2
SQ_eff = SQ_Pmax;

% Export the SQ data
writematrix([SQ_Voltage SQ_Current_Density], strcat(Output_Data_Directory, ...
    "Shockley-Queisser Model JV Curve.txt"), 'Delimiter', 'tab');
writematrix([SQ_Voc; SQ_Jsc; SQ_FF; SQ_eff], strcat(Output_Data_Directory, ...
    "Shockley-Queisser Parameters.txt"), 'Delimiter', 'tab');
%% Section 2: Calculate efficiency limits with extended Tiedje-Yablonovitch
% (TY) Model (includes Auger recombination) and defect-assisted
% Shockley-Reed-Hall (SRH) recombination

% Data output tables
% Parameters
Voc = zeros(length(L_Range), length(SRH_tau_range));
Jsc = zeros(length(L_Range), length(SRH_tau_range));
FF = zeros(length(L_Range), length(SRH_tau_range));
Vmax = zeros(length(L_Range), length(SRH_tau_range));
Imax = zeros(length(L_Range), length(SRH_tau_range));
eff = zeros(length(L_Range), length(SRH_tau_range));

% Recombination components
AugerComponent = zeros(length(L_Range), length(SRH_tau_range));
FreeCarrierComponent = zeros(length(L_Range), length(SRH_tau_range));
InternalLuminescenceComponent = zeros(length(L_Range), length(SRH_tau_range));
ReabsorptionComponent = zeros(length(L_Range), length(SRH_tau_range));
ExternalEmissionComponent = zeros(length(L_Range), length(SRH_tau_range));
SRHComponent = zeros(length(L_Range), length(SRH_tau_range));

% Lifetimes
Auger_lifetime = zeros(length(L_Range), length(SRH_tau_range));
Radiative_lifetime = zeros(length(L_Range), length(SRH_tau_range));
SRH_lifetime = zeros(length(L_Range), length(SRH_tau_range));

% Supplementary Note 1, Equation 3
Black_Body_Photon_Flux = (2 * Material_n2) / (speed_of_light ^ 2 * ...
    planck_eV ^ 3) * (New_Energy_eV_range .^ 2) ./ (exp(New_Energy_eV_range ...
    ./ ((Boltzmann_constant / electron_charge) * T)) - 1);

thickness_interest_index = 0;
for thickness_index = 1:length(L_Range)

    % Print current thickness and time
    L = L_Range(thickness_index);
    fprintf(['Current Thickness: ' num2str(L), ' nm (Started), ']);
    Time1 = datetime('now');
    fprintf(['Current Time Now: ' datestr(Time1) '.\n']);

    if ismember(L, L_Interest)
        [~, ~, ~] = mkdir(strcat("L = ", num2str(L), " nm Data"));
        L_Interest_Directory = strcat("L = ", num2str(L), " nm Data/");
    end

    % Supplementary Note 1, Equation 2
    Absorbance = Material_Alpha2 ./ (Material_Alpha2 + Material_Alpha1 + (1 ./ (Material_4n2 .* L)));
    
    f_range = [0; linspace(0.8, 1, 100).'];
    if ismember(L, L_Interest) % use whole range for thicknesses of interest
        f_range = linspace(0, 1, 1000).';
    end
    
    Current_Density = zeros(length(f_range), length(SRH_tau_range));
    Voltage = zeros(length(f_range), length(SRH_tau_range));
    
    for tau_index = 1:length(SRH_tau_range)
        SRH_tau = SRH_tau_range(tau_index);

        % Supplementary Note 1, Equation 4
        Jsc(thickness_index, tau_index) = 10 ^ 3 * 10 ^ -4 * electron_charge *...
            trapz(New_Wavelength_Range, Spectrum .* Absorbance ./ New_Energy_J_Range);
        % integrand units: (J / s) / (nm * m ^ 2) / J = 1 / (nm * m ^ 2 * s)
        % units: (photon * C) / (nm * m ^ 2 * s) = A / m ^ 2 = 10 ^ 3 * 10 ^ -4 * mA / cm ^ 2

        for fraction_index = 1:length(f_range)
            f = f_range(fraction_index);

            syms x;
            % Supplementary Note 1, Equation 11
            % units (first term): (nm ^ -1 + nm ^ -1) * 10 ^ 7 nm / cm = 1 / cm
            % units (second term): (10 ^ -4 m ^ 2 / cm ^ 2) * (1 / (s * m ^ 2 * eV)) * eV = cm ^ -2 / s
            % units (third term): cm ^ 6 * s ^ -1 * cm ^ -9 = cm ^ -3 / s
            % units (right hand side): ((C / s) / cm ^ 2) / (C * nm * 10 ^ -7 cm / nm) = cm ^ -3 / s
            if SRH_tau == Inf % No SRH
                eqn_RadiativeAugerFreeSRH = 4 * pi * (10 ^ 7 * (Material_Alpha1 * ...
                    exp(x / 2) + 1 / (Material_4n2 .* L))) * exp(x) * 10 ^ -4 * ...
                    trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* ...
                    Absorbance)) + Material_C * Material_ni ^ 3 * exp((3 / 2) * ...
                    x) == (Jsc(thickness_index, tau_index) * 10 ^ -3) / ...
                    (electron_charge * L * 10 ^ -7) * (1 - f);
            else
                eqn_RadiativeAugerFreeSRH = 4 * pi * (10 ^ 7 * (Material_Alpha1 * ...
                    exp(x / 2) + 1 / (Material_4n2 .* L))) * exp(x) * 10 ^ -4 * ...
                    trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* ...
                    Absorbance)) + Material_C * Material_ni ^ 3 * exp((3 / 2) * ...
                    x) + Material_ni * exp (x / 2) / SRH_tau == ...
                    (Jsc(thickness_index, tau_index) * 10 ^ -3) / (electron_charge * ...
                    L * 10 ^ -7) * (1 - f);
            end

            X_RadiativeAugerFreeSRH = vpasolve(eqn_RadiativeAugerFreeSRH, x);

            Current_Density(fraction_index, tau_index) = f * Jsc(thickness_index, tau_index);
            if f == 1
                if any(Voltage(:, tau_index) < 0)
                    Voltage(fraction_index, tau_index) = Voltage(fraction_index - 1, tau_index);
                    Jsc(thickness_index, tau_index) = interp1(Voltage(1:end - 1, tau_index), ...
                        Current_Density(1:end - 1, tau_index), 0);
                else
                    Voltage(fraction_index, tau_index) = 0;
                end
            else
                Voltage(fraction_index, tau_index) = (Boltzmann_constant * ...
                    T * X_RadiativeAugerFreeSRH) / (electron_charge);
            end

            if f == 0
                Voc(thickness_index, tau_index) = (Boltzmann_constant * T * ...
                    X_RadiativeAugerFreeSRH) / (electron_charge);
            end
        end

        % Update the table with parameters
        FF(thickness_index, tau_index) = max(Current_Density(:, tau_index) .* Voltage(:, tau_index)) /...
            (Voc(thickness_index, tau_index) * Jsc(thickness_index, tau_index));
        for max_index = 1:length(Voltage(:, tau_index))
           if Voltage(max_index, tau_index) * Current_Density(max_index, tau_index) == max(Voltage(:, tau_index) .* Current_Density(:, tau_index))
               Vmax(thickness_index, tau_index) = Voltage(max_index, tau_index);
               Imax(thickness_index, tau_index) = Current_Density(max_index, tau_index);
           end
        end
        eff(thickness_index, tau_index) = max(Current_Density(:, tau_index) .* Voltage(:, tau_index)) / (trapz(New_Wavelength_Range, Spectrum) * 10 ^ 3 * 10 ^ -4) * 100;

        % units: mA / cm ^ 2
        AugerComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * Material_C * Material_ni ^ 3 * ...
            exp((3 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)));
        FreeCarrierComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
            exp((3 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T))) * ...
            (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 * Material_Alpha1)));
        InternalLuminescenceComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
            exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
            (Material_Alpha1 * exp((1 / 2) * (electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
            trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7)) + ...
            trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)) + ...
            (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance)));
        ReabsorptionComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
            exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
            (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)));
        ExternalEmissionComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * 4 * pi * 10 ^ -4 * ...
            exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)) * ...
            (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance));
        SRHComponent(thickness_index, tau_index) = electron_charge * L * 10 ^ 3 * 10 ^ -7 * Material_ni * ...
            exp ((1 / 2) * ((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T))) / SRH_tau;
       
        % units: s
        Auger_lifetime(thickness_index, tau_index) = 1 / (Material_C * Material_ni ^ 2 * exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)));
        Radiative_lifetime(thickness_index, tau_index) = 1 ./ ((Material_ni * sqrt(exp((electron_charge * Vmax(thickness_index, tau_index)) / (Boltzmann_constant * T)))) .* ...
            (((4 * pi * 10 ^ -4 * (trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 * Material_Alpha1)) + ...
            trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance * 10 ^ 7 .* Material_Alpha2)) + ...
            (1 / (Material_4n2 .* L * 10 ^ -7)) * trapz(flip(New_Energy_eV_range), flip(Black_Body_Photon_Flux .* Absorbance))))) ./ ...
            (Material_ni ^ 2)));
        SRH_lifetime(thickness_index, tau_index) = SRH_tau;
    end
    
    % Generate JV Curve (Supplementary Figure 3)
    if ismember(L, L_Interest)
        SQ_JV = readmatrix(strcat(Output_Data_Directory, ...
            "Shockley-Queisser Model JV Curve.txt"));
        SQ_Parameters = readmatrix(strcat(Output_Data_Directory, ...
            "Shockley-Queisser Parameters.txt"));

        FIGURE_COUNTER = FIGURE_COUNTER + 1;
        figure(FIGURE_COUNTER);

        hold on;
        for j = 1:length(SRH_tau_range)
            plot(Voltage(:, j), Current_Density(:, j), '-', 'color', Colors(j, :), 'linewidth', 3);
        end
        plot(SQ_JV(:, 1), -SQ_JV(:, 2), '-', 'color', '[0 0 0]', 'linewidth', 3);

        xlabel('Voltage (V)');
        ylabel('Current density (mA cm^{-2})');
        hTitle = title(strcat("JV Curve, Thickness = ", num2str(L), " nm"));
        set(gca, 'FontSize', 20);
        set(hTitle, 'FontSize', 32);
        set(gca, 'linewidth', 2);
        xlim([-0.1 * SQ_Parameters(1) 1.25 * SQ_Parameters(1)]);
        ylim([-0.1 * SQ_Parameters(2) 1.25 * SQ_Parameters(2)]);
        set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);
        
        % Export the plot and the data
        savefig(strcat(L_Interest_Directory, "JV Curve, ", num2str(L), " nm.fig"));
        saveas(gcf, strcat(L_Interest_Directory, "JV Curve, ", num2str(L), " nm.jpg"));
        close;
        writematrix([Voltage Current_Density], strcat(L_Interest_Directory, ...
            "JV Data, ", num2str(L), " nm.txt"), 'Delimiter', 'tab');
    end

    % Generate spectral dependence of the luminescent emission rate (Figure 4)
    if ismember(L, L_Interest)
        FIGURE_COUNTER = FIGURE_COUNTER + 1;
        figure(FIGURE_COUNTER);

        % Sum of all 3 terms (Total Internal)
        plot(New_Energy_eV_range, 10 ^ -6 * (4 * pi * 10 ^ -4 * Black_Body_Photon_Flux .* ...
            Absorbance * 10 ^ 7 * Material_Alpha1 + 4 * pi * 10 ^ -4 * Black_Body_Photon_Flux .* ...
            Absorbance * 10 ^ 7 .* Material_Alpha2 + (1 / (Material_4n2 * L * 10 ^ -7)) * ...
            4 * pi * 10 ^ -4 * Black_Body_Photon_Flux .* Absorbance), ...
            '-', 'color', '[0.25 0.25 0.25]', 'linewidth', 3);
        hold on;

        % 2nd term (Reabsorbtion)
        plot(New_Energy_eV_range, 10 ^ -6 * (4 * pi * 10 ^ -4 * Black_Body_Photon_Flux .*...
            Absorbance * 10 ^ 7 .* Material_Alpha2), '-', 'color', '[0.00 0.45 0.74]', 'linewidth', 3);
        
        % 3rd term (External Emission)
        plot(New_Energy_eV_range, 10 ^ -6 * ((1 / (Material_4n2 * L * 10 ^ -7)) * 4 * pi *...
            10 ^ -4 * Black_Body_Photon_Flux .* Absorbance), '-', 'color', '[0.85 0.33 0.10]', 'linewidth', 3);

        xlabel('Energy (eV)');
        ylabel({'Luminescent emission'; '(10^6 photons cm^{-3} s^{-1} eV^{-1})'});
        hTitle = title(strcat("Luminescent Emission Rate, Thickness = ", num2str(L), " nm"));
        set(gca, 'FontSize', 20);
        set(hTitle, 'FontSize', 32);
        set(gca, 'linewidth', 2);
        xlim([1.0 1.6]);
        set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

        savefig(strcat(L_Interest_Directory, "Luminescent Emission.fig"));
        saveas(gcf, strcat(L_Interest_Directory, "Luminescent Emission.jpg"));
        close;
    end

    % Print current thickness and time
    fprintf(['Current Thickness: ' num2str(L), ' nm (Finished), ']);
    Time2 = datetime('now');
    fprintf(['Current Time Now: ' datestr(Time2) '.\n']);

    if ismember(L, L_Interest)
        Time_Logger_Interest = Time2 - Time1;

        % Move to the next thickness of interest
        thickness_interest_index = thickness_interest_index + 1;
    else
        Time_Logger_Regular = Time2 - Time1;
    end

    % Print number of thicknesses remaining
    remaining_lengths = length(L_Range) - thickness_index;
    remaining_lengths_interest = length(L_Interest) - thickness_interest_index;
    fprintf(['Thicknesses Remaining: ' num2str(remaining_lengths)]);
    
    % Print approximate time remaining
    if thickness_interest_index > 0
        % More accurate prediction based on time it takes to calculate
        % thicknesses of interests
        remaining_time = (remaining_lengths - remaining_lengths_interest) * ...
            Time_Logger_Regular + remaining_lengths_interest * Time_Logger_Interest;
    else
        remaining_time = remaining_lengths * Time_Logger_Regular;
    end
    [hours, minutes, seconds] = hms(remaining_time);
    fprintf([', ' 'Approximate Time Remaining: ' num2str(hours) ' Hours, ' ...
        num2str(minutes) ' Minutes, ' num2str(seconds) ' Seconds\n']);
end

% Export the data
writematrix([L_Range Voc], strcat(Output_Data_Directory, "Voc.txt"), 'Delimiter', 'tab');
writematrix([L_Range Jsc], strcat(Output_Data_Directory, "Jsc.txt"), 'Delimiter', 'tab');
writematrix([L_Range FF], strcat(Output_Data_Directory, "FF.txt"), 'Delimiter', 'tab');
writematrix([L_Range Vmax], strcat(Output_Data_Directory, "Vmax.txt"), 'Delimiter', 'tab');
writematrix([L_Range Imax], strcat(Output_Data_Directory, "Imax.txt"), 'Delimiter', 'tab');
writematrix([L_Range eff], strcat(Output_Data_Directory, "eff.txt"), 'Delimiter', 'tab');

writematrix([L_Range AugerComponent], strcat(Output_Data_Directory, "Component_Auger.txt"), 'Delimiter', 'tab');
writematrix([L_Range FreeCarrierComponent], strcat(Output_Data_Directory, "Component_FreeCarrier.txt"), 'Delimiter', 'tab');
writematrix([L_Range InternalLuminescenceComponent], strcat(Output_Data_Directory, "Component_InternalLuminescence.txt"), 'Delimiter', 'tab');
writematrix([L_Range ReabsorptionComponent], strcat(Output_Data_Directory, "Component_Reabsorption.txt"), 'Delimiter', 'tab');
writematrix([L_Range ExternalEmissionComponent], strcat(Output_Data_Directory, "Component_ExternalEmission.txt"), 'Delimiter', 'tab');
writematrix([L_Range SRHComponent], strcat(Output_Data_Directory, "Component_SRH.txt"), 'Delimiter', 'tab');

writematrix([L_Range Auger_lifetime], strcat(Output_Data_Directory, "Lifetime_Auger.txt"), 'Delimiter', 'tab');
writematrix([L_Range Radiative_lifetime], strcat(Output_Data_Directory, "Lifetime_Radiative.txt"), 'Delimiter', 'tab');
writematrix([L_Range SRH_lifetime], strcat(Output_Data_Directory, "Lifetime_SRH.txt"), 'Delimiter', 'tab');
%% Section 3: Plotting the data (Open Circuit Voltage, Short Circuit
% Current Density, Fill Factor, and Power Conversion Efficiency)

SQ_Parameters = readmatrix(strcat(Output_Data_Directory, ...
    "Shockley-Queisser Parameters.txt"));

Text_Filename = ["Voc"; "Jsc"; "FF"; "eff"];

Y_Labels = ["Open-circuit voltage (V)"; "Short circuit current density (mA cm^{-2})"; ...
    "Fill factor (%)"; "Power conversion efficiency (%)"];

Titles = ["Open Circuit Voltage"; "Short Circuit Current Density"; ...
    "Fill Factor"; "Power Conversion Efficiency"];

for i = 1:length(Text_Filename)
    Parameters = readmatrix(strcat(Output_Data_Directory, Text_Filename(i), ".txt"));
    FIGURE_COUNTER = FIGURE_COUNTER + 1;
    figure(FIGURE_COUNTER);

    hold on;
    for j = 1:length(SRH_tau_range)
        semilogx(Parameters(:, 1), Parameters(:, j + 1), '-', 'color', Colors(j, :), 'linewidth', 3);
    end
    semilogx([Parameters(1, 1) Parameters(end, 1)], [SQ_Parameters(i) SQ_Parameters(i)], '-', 'color', '[0 0 0]', 'linewidth', 3);

    xlabel('Thickness (nm)');
    ylabel(Y_Labels(i));
    hTitle = title(strcat(Material_Name, " ", Titles(i), " vs. Thickness"));
    set(gca, 'FontSize', 20);
    set(hTitle, 'FontSize', 32);
    set(gca, 'XScale', 'log');
    set(gca, 'linewidth', 2);
    xlim([5 1000]);
    set(gcf, 'Units', 'normalized', 'Position', [0 0 1 1]);

    savefig(strcat(Figures_fig_Directory, Titles(i), " vs. Thickness.fig"));
    saveas(gcf, strcat(Figures_jpg_Directory, Titles(i), " vs. Thickness.jpg"));
    close;
end