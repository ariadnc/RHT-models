
%****************************************************************************************
%
%               MATLAB CODE WRITEN BY ARIADNA CRESPO AND CHATI (XD)  
%
% Discription: Extract info of th ename of the file sto anlysis
%=========================================================================================

function [fuel, pressure_val, O2_pct, dil_pct, diluent] = parseFlameFilename(filename)
    [~, name, ~] = fileparts(filename);  % strip extension
    pattern = '([A-Za-z0-9]+)_([0-9\.]+)([a-zA-Z]+)_([0-9]+)_O2_([0-9]+)_([A-Za-z0-9]+)';
    tokens = regexp(name, pattern, 'tokens');

    if isempty(tokens)
        error("Filename does not match expected pattern: %s", filename);
    end

    parts = tokens{1};

    fuel         = parts{1};
    pressure_str = [parts{2} parts{3}];  % e.g., "1atm"
    pressure_val = str2double(parts{2});  % e.g., 1
    O2_pct       = str2double(parts{4});
    dil_pct      = str2double(parts{5});
    diluent      = erase(parts{6}, '.csv');  % extra safety (probably not needed anymore)
end
