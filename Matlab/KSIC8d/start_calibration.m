clear;
clc;
rng('default');    % random number generator

%% Assigned Parameters


S  = 300;         % 1/2 number of sectors
L  = 1;            % inelastic labor supply in Home
Ls = 1;            % inelastic labor supply in Foreign

N  = 150;          % maximum number of producers per sector

gamma = 10;                   % within-sector elasticity of substitution   
theta = 1.25;                 % across-sector elasticity of substitution
zeta = 0.02;


%% Calibration Results - Bertrand Model

results_bertrand = [];
results_cournot = [];
comp_trade_elasticity = 1;

for xi_x = 5
    for xi_z = 5
        for FD = 0.002
            for FX = 0.01
                for tau = 0.13
                    for rho = 1
                          
                        par_value = [xi_x, xi_z, zeta, FD, FX, tau, rho];

                        equilibrium_logn_bertrand;
                        model_moments_calibration;
                        results_bertrand = [results_bertrand; [par_value,all_model_moments,sigma]];

%                         equilibrium_cournot;
%                         model_moments_calibration;
%                         results_cournot = [results_cournot; [parameters, all_model_moments']];

                    end
                end
            end
        end
    end
end

%xlswrite('calibration_results_bertrand.xlsx', results_bertrand);
%xlswrite('calibration_results_cournot.xlsx', results_cournot);

