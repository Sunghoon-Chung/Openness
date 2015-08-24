
%% Productivity backout

rng(0);   % reset random number generator

ws = tds/sum(tds);          % share of sector s in the total demand (in home country)
mud = mu;
muf = ((gamma-1)/gamma - (1/theta - 1/gamma).*wfs).^(-1);

if tariff==1
    tar = ot;                 % tariff rate
else
    tar = zeros(S,N);
end


%% initial guess

ads = lp;            % labor productivity as initial value of ads 
afs = ones(S,1);


%% iterate until convergence

for i = 1:500

    ads_new = mud.*wds.^(1/(gamma-1))./(ws.*mean((mud.*wds.^(1/(gamma-1))./ads).^(1-theta))).^(1/(1-theta));
    afs_new = (1./(1-tar).*muf.*wfs.^(1/(gamma-1)))./(ws.*mean(muf.*wfs.^(1/(gamma-1))./(afs.*(1-tar))).^(1-theta)).^(1/(1-theta));
    error = norm([ads_new-ads, afs_new-afs]);        %% norm() = Euclidean norm
    %fprintf('%4i %6.2e \n',[i, error]);

    if error < 1e-7, break, end
    
    ads = real(ads_new);
    afs = real(afs_new);
    
end

ads = ads./min(ads);
scatter(ads, lp);
afs = afs./min(afs);
histogram(afs);
scatter(ads, afs);
[rho, pval] = corr(ads, afs, 'type', 'Kendall');


%% Method 2

ksdensity(ads);
xlabel('Domestic Productivity');
ads_param = mle(ads, 'distribution', 'Lognormal');
afs = lognrnd(ads_param(1), ads_param(2), S,1);
ksdensity(afs);
xlabel('Foreign Productivity');

[adss, ads_id] = sort(ads, 1, 'descend');
[afss, afs_id] = sort(afs, 1, 'descend');
scatter(ads, afs);
[rho, pval] = corr(adss, afss, 'type', 'Spearman');


