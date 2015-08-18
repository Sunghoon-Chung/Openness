clear;
rng('default');


% Assigned Parameters and raw data

year = 2005; 

io_data = strcat('D:\Copy\Openness\Analysis\Matlab\io_data_',num2str(year),'.mat');
load(io_data)

b1 = -0.1092;       % slope coefficient and intercept in regression of inverse markups on market shares
b0 = 0.6838;

S = n_sector(1);

%gamma = 1/(1-b0);                                    % within-sector elasticity of substitution   
gamma = 8.335;
theta = (1/gamma - b1/b0*(1-1/gamma))^(-1);     % across-sector elasticity of substitution

L  = 1;            % inelastic labor supply in Home
Ls = 1;            % inelastic labor supply in Foreign

N  = 1;            % fixed N, no geometric distribution of n(s)

tariff = 1;        % 1 if yes, 0 if no



% Calibrated Parameters

FD    = 0;                                  % fixed cost of domestic operations
FX    = 0;                                  % fixed cost of export operations

tau   = 0.13;                                  % net trade cost 


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
     
scatter(omegaH(1:S), muH(1:S));
scatter(wds, mud);        

