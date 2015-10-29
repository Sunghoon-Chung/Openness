clear;
rng('default');


% Assigned Parameters and raw data

year=1973;

timevarying = 0;
if timevarying==1
    if     year==1970,  b1 = -0.1225;   b0 = 0.6291;   % slope coefficient and intercept in regression of inverse markups on market shares
    elseif year==1973,  b1 = -0.1613;   b0 = 0.6693;
    elseif year==1980,  b1 = -0.1361;   b0 = 0.7235;
    elseif year==1990,  b1 = -0.1454;   b0 = 0.7482;
    elseif year==2003,  b1 = -0.2026;   b0 = 0.7670;
    elseif year==2005,  b1 = -0.1889;   b0 = 0.7922;
    end
else b1= -0.1454; b0= 0.7482;
end
    
    
io_data = strcat('D:\Copy\Openness\Analysis\Matlab\IOtable\io_data_',num2str(year),'.mat');
load(io_data)

S = n_sector(1);

%gamma = 1/(1-b0);                                    % within-sector elasticity of substitution   
gamma = 10;
theta = (1/gamma - b1/b0*(1-1/gamma))^(-1);     % across-sector elasticity of substitution

L  = 1;            % inelastic labor supply in Home
Ls = 1;            % inelastic labor supply in Foreign

N  = 1;            % fixed N, no geometric distribution of n(s)

tariff = 1;        % 1 if yes, 0 if no



% Calibrated Parameters

FD    = 0;                                  % fixed cost of domestic operations
FX    = 0;                                  % fixed cost of export operations

tau   = 0.12;                                  % net trade cost 


tausave  = tau;   % save since need later
FXsave   = FX; 

% Uncomment the following to report more statistics

%model_moments;
%markup_moments;
     
%break 
 


% Cumputing Equilibrium with Trade

fprintf('\n');
display('*** Computing Equilibrium with Trade ***')
fprintf('\n');

tau = tausave;
FX  = FXsave;
comp_trade_elasticity = 1;
        
equilibrium;
        
 
% fprintf('A increase (gains from trade), *100      = %7.3f \n', log(A/Asave)*100);  
% fprintf('Due to markups                           = %7.3f \n', Alosssave - log(Aeff/A)*100);  
% fprintf('Trade elasticity                         = %7.3f \n', sigma);  
% fprintf('Import share                             = %7.3f \n', impshare);
% fprintf('Fraction exporters                       = %7.3f \n', agg_fexporters);  

% Uncomment the following to report more statistics

%model_moments;
%markup_moments;
%domestic_and_import_markups;         

