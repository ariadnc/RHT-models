%=========================================================================================
% Created on Tue Aug 26 16:09:24 2025
%
% @author: ariad
%=========================================================================================

folders = {
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\CH4_1atm\CO2', ...
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\CH4_1atm\N2', ...
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\H2_1atm\CO2', ...
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\H2_1atm\N2', ...
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\H2_10atm\CO2', ...
    'C:\Users\ariad\OneDrive\Documentos\MATLAB\H2_10atm\N2'...
};

% === Constants ===
L = 0.01;                  % Path length (m)
sigma = 5.670374419e-11;   % kW/m^2·K^4
T_amb = 700;               % K
C_part = 1186;             % soot coefficient 

% === Loop over each folder ===
for folder_idx = 1:length(folders)
    folder_path = folders{folder_idx};
    files = dir(fullfile(folder_path, '*.csv'));

    % Containers for each folder
    all_eps      = {};
    all_kabs     = {};
    all_x        = {};
    O2_values    = [];
    diluent_name = "";
    all_k_total  = {};
    all_Q_rad_vol= {};
    all_qr_OTA   = {};   

    % === Loop through flame files in folder ===
    for f = 1:length(files)
        filename = files(f).name;
        filepath = fullfile(folder_path, filename);

        data = readtable(filepath);

        x_cm = data.("Distance_cm_");
        T = data.("Temperature_K_");
        x_h2o = data.("Mole_fraction_H2O__");

        if ismember("Mole_fraction_CO2__", data.Properties.VariableNames)
            x_co2 = data.("Mole_fraction_CO2__");
        else
            x_co2 = ones(size(x_h2o)) * 1e-10; % So no error when dividing /0
        end

        try
            [fuel, pressure, O2_pct, ~, diluent] = parseFlameFilename(filename);
            O2_values(end+1) = O2_pct;
            diluent = erase(diluent, 'csv');
            diluent_name = diluent;  % same for all in this folder
        catch
            warning("Skipping file due to parsing issue: %s", filename);
            continue;
        end

        % Fuel properties 
        if strcmp(fuel, "CH4")
            LHV = 50e3;         % kJ/kg
            rho_fuel = 1.1;     % kg/m^3
            v = 25/100;         % m/s
        elseif strcmp(fuel, "H2")
            LHV = 120e3;        % kJ/kg
            rho_fuel = 0.085;   % kg/m^3
            v = 137.5/100;      % m/s
        end

        % === Soot column (volume fraction) ===
        if ismember("Particle_volume_fraction_cm3_cm3_", data.Properties.VariableNames)
            fv = data.("Particle_volume_fraction_cm3_cm3_");
        else
            fv = zeros(size(T));
        end

        % === Compute emissivity and absorption ===
        P_total = pressure;   % atm
        n = length(T);
        eps_arr     = zeros(n, 1);
        k_abs_arr   = zeros(n, 1);
        k_soot_arr  = zeros(n, 1);
        k_total_arr = zeros(n, 1);
        Q_rad_vol_arr = zeros(n, 1);

        for i = 1:n
            T_i   = T(i);
            xh2o  = x_h2o(i);
            xco2  = max(x_co2(i), 1e-6);
            Mr    = xh2o / xco2;

            if Mr > 4,   Mr = 4;   xco2 = xh2o / Mr; end
            if Mr < 0.1, Mr = 0.1; xh2o = xco2 * Mr; end
            if isnan(T_i), eps_arr(i) = NaN; k_abs_arr(i) = NaN; continue; end

            [~, ~, eps] = Bordbar_PWSGG_MatlabFunction_IJHMT2021(P_total, T_i, xh2o, xco2, L);
            eps = min(max(eps, 0), 1);
            eps_arr(i) = eps;

            k_abs = -log(max(1 - eps, 1e-10)) / L;
            k_abs_arr(i) = k_abs;
            
            k_soot = C_part * fv(i) * T(i);
            k_soot_arr(i) = k_soot;
            
            k_total = k_abs + k_soot;
            k_total_arr(i) = k_total;

            Q_rad_vol_arr(i) = 4 * k_total * sigma * (T_i^4 - T_amb^4);
            
        end

        % OTA source and flux (integrate along x)
        dx_cm = mean(diff(x_cm)); 
        Sr_OTA = Q_rad_vol_arr;
        qr_OTA = cumtrapz(x_cm/100, Sr_OTA);  % kW/m^2

        % Chem heat input model
        m_dot  = rho_fuel * v;     % kg/m^2/s
        q_chem = m_dot * LHV;      % kW/m^2

        Xrad = qr_OTA(end) / q_chem;

        OTA_Sr_data(key) = Sr_OTA;
        OTA_qr_data(key) = qr_OTA;

        all_eps{end+1} = eps_arr(:);
        all_kabs{end+1} = k_abs_arr(:);
        all_x{end+1} = x_cm(:);
        all_k_total{end+1} = k_total_arr(:);
        all_Q_rad_vol{end+1} = Q_rad_vol_arr(:);
        all_qr_OTA{end+1} = qr_OTA(:);
    end
    
end

