

fprintf('...calculating model moments... \n');

%%%%% A. Domestic market share moments   

omegad = omegaH;
omegad = omegaH./repmat(sum(omegaH),N,1);
omegad = omegad(:,sum(omegaH)>0); 

[NN,SS] = size(omegad);

omegad(omegad==0) = NaN; %% mimic data structure (use of NaN instead of zeros)

sharedom = omegad;       %% domestic market shares

hh            = nansum(sharedom.^2);
highest_share = nanmax(sharedom);

share         = sharedom(:); share = share(share>0);


for s = 1:SS,
    
    share90(s)          = prctile(sharedom(:,s),90);
    share75(s)          = prctile(sharedom(:,s),75);
    share50(s)          = prctile(sharedom(:,s),50);
    share25(s)          = prctile(sharedom(:,s),25);
    share10(s)          = prctile(sharedom(:,s),10);
        
    sdshare(s)          = nanstd(sharedom(:,s));
        
    inversehh(s)        = 1./nansum(sharedom(:,s).^2);
        
    topshare(s)         = nanmax(sharedom(:,s));
        
    ndom(s)             = sum((sharedom(:,s)>0));
                
    ihh(s)              = 1./nansum(sharedom(:,s).^2);
              
end

%%% collate

HH            = [hh];
Highest_Share = [highest_share];
Share         = [share];
    
Share90       = [share90];
Share75       = [share75];
Share50       = [share50];
Share25       = [share25];
Share10       = [share10];
    
Sdshare       = [sdshare];
    
Ndom          = [ndom];
    
IndustryIHH   = [inversehh];
IndustryTop   = [topshare];


%%%%% B. Size distribution moments

sales = pH.*yH; %% domestic only
labor = lH;

sales_s = sum(sales)';
labor_s = sum(labor)';

labor = labor(sales>0);
sales = sales(sales>0);

data = sortrows([sales, labor],1);

fRtopY1  = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); %% size distribution of establishments
fRtopY5  = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); %% (data counterpart is value-added)
fRtopL1  = sum(data(end-floor(length(data)*0.01)+1:end,2))/sum(data(:,2));
fRtopL5  = sum(data(end-floor(length(data)*0.05)+1:end,2))/sum(data(:,2));

data  = sortrows([sales_s, labor_s], 1);

fRtopR1_s = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); %% size distribution of sectors
fRtopR5_s = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); %% (data counterpart is sales)
fRtopL1_s = sum(data(end-floor(length(data)*0.01)+1:end,2))/sum(data(:,2));
fRtopL5_s = sum(data(end-floor(length(data)*0.05)+1:end,2))/sum(data(:,2));


%%%%% C. Export/Import moments

fexporters     = sum(omegaH>0 & omegaF>0)./sum(omegaH>0);
agg_fexporters = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0); 

imports_s      = sum(pF.*yF);  
import_shares  = imports_s./sum(pF.*yF+pH.*yH); 
agg_impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:));

impshare   = sum(omegaF)';
expshare   = sum(omegaHs)';
ttrade     = impshare + expshare;


%%%%% Regression that identifies theta from inv markups and market shares

tempX = [ones(size(omegaH(:))),omegaH(:)];
tempY = 1./muH(:);

temp_beta  = inv(tempX'*tempX)*tempX'*tempY;

beta_invmarkup_marketshare = temp_beta(2);



%%%%% GL index

intraindex = 1 - mean(sj(ttrade>0)'.*abs(impshare(ttrade>0)-expshare(ttrade>0))./(impshare(ttrade>0) + expshare(ttrade>0)))./mean(sj(ttrade>0)');

%%%%% correlations/regressions of sector import, export shares and sales shares

sector_imports       = sum(pF.*yF);
total_imports        = sum(sector_imports);
sector_share_imports = sector_imports/total_imports;

sector_exports       = sum(pHs.*yHs);
total_exports        = sum(sector_exports);
sector_share_exports = sector_exports/total_exports;

sector_sales         = sum(pH.*yH); 
total_sales          = sum(sector_sales);
sector_share_sales   = sector_sales/total_sales;


temp = corrcoef([sector_share_imports',sector_share_exports',sector_share_sales']);

corr_share_exports_share_sales = temp(2,3); %% correlation of sector share exports with sector share of sales
corr_share_imports_share_sales = temp(1,3); %% correlation of sector share imports with sector share of sales

tempX = [ones(size(sector_share_imports))',sector_share_sales'];
tempY = [sector_share_imports'];

temp_beta = inv(tempX'*tempX)*tempX'*tempY;

beta_share_imports_share_sales = temp_beta(2); %% regression coefficient of sector share imports on sector share *sales*

tempX = [ones(size(sector_share_imports))',sector_share_exports'];
tempY = [sector_share_imports'];

temp_beta = inv(tempX'*tempX)*tempX'*tempY;

beta_share_imports_share_exports = temp_beta(2); %% regression coefficient of sector share imports on sector share *exports*

%%%%% Import share and concentration regressions

% New regressions of import shares on HH 
 
tempX = [ones(SS,1), HH'];
tempY = [sector_share_imports(sum(omegaH)>0)'];
 
temp_betaHH = inv(tempX'*tempX)*tempX'*tempY;
 
% New regressions of import share on top share: 
 
tempX = [ones(SS,1), topshare'];
tempY = [sector_share_imports(sum(omegaH)>0)'];
 
temp_betaTOP = inv(tempX'*tempX)*tempX'*tempY;

%%%% Regression to identify correlation in x_i

ind_sales=sector_sales(sum(omegaH)>0)'+sector_imports(sum(omegaH)>0)';
import_pen=sector_imports(sum(omegaH)>0)'./ind_sales;
tempY = [import_pen];
tempX = [ones(SS,1), HH'];

temp_beta= inv(tempX'*tempX)*tempX'*tempY;

beta_imp_HH = temp_beta(2);

% Blocks are to make it easier to move them around if required

% Block #1

meaninvHH         = nanmean(1./HH); 
medianinvHH       = nanmedian(1./HH);
meanmaxshare      = nanmean(Highest_Share);
medianmaxshare    = nanmedian(Highest_Share);          

Block1            = [meaninvHH,medianinvHH,meanmaxshare,medianmaxshare];

% Block #2 

meanshare         = nanmean(Share);
sdshare           = nanstd(Share); 
medianshare       = nanmedian(Share);
ppshare           = prctile(Share,[75,95,99]);         

Block2            = [meanshare,sdshare,medianshare,ppshare];

% Block #3


Block3            = [fRtopR1_s,fRtopR5_s,fRtopL1_s,fRtopL5_s];

% Block #4 

ppndom            = prctile(Ndom,[10,25,50,75,90,95]); 

Block4            = [ppndom];

% Block #5

Block5            = [fRtopY1,fRtopY5,fRtopL1,fRtopL5];

% Block #6

Block6            = [agg_fexporters,agg_impshare,beta_invmarkup_marketshare,aweight,intraindex,corr_share_exports_share_sales,corr_share_imports_share_sales,beta_share_imports_share_sales,beta_share_imports_share_exports];

% Block #7

mean90to10        = nanmean(Share90-Share10);
median90to10      = nanmedian(Share90-Share10);
mean75to25        = nanmean(Share75-Share25);
median75to25      = nanmedian(Share75-Share25);
meanSDshare       = nanmean(Sdshare);
medianSDshare     = nanmedian(Sdshare);                

Block7            = [mean90to10,median90to10,mean75to25,median75to25,meanSDshare,medianSDshare];

% Block #8

ppindustryihh     = prctile(IndustryIHH,[10,25,50,75,90,95]); 

Block8            = [ppindustryihh];

% Block #9

ppindustrytop     = prctile(IndustryTop,[10,25,50,75,90,95]); 

Block9            = [ppindustrytop];


all_model_moments = [Block1,Block2,Block3,Block4,Block5,Block6,Block7,Block8,Block9]';



%%%%% load data moments from saved file

load saved_data_moments_march2015.mat; % vector of data moments, ordered same as model moments above
                                     % particular care should be taken if there is any
                                     % reordering of moments or change in
                                     % moment selection
                                  

if length(all_model_moments)~=length(all_data_moments);
    
    display('problem with number of model and data moments!')
    break
end


        

%%%%% list of moments

display('                                         Model    Data')
fprintf('\n');
fprintf('mean HH inverse                      = %7.3f %7.3f \n',nanmean(1./HH)             ,all_data_moments(1));  
fprintf('median HH inverse                    = %7.3f %7.3f \n',nanmedian(1./HH)           ,all_data_moments(2));
fprintf('mean highest share                   = %7.3f %7.3f \n',nanmean(Highest_Share)     ,all_data_moments(3));  
fprintf('median highest share                 = %7.3f %7.3f \n',nanmedian(Highest_Share)   ,all_data_moments(4));
fprintf('\n');
fprintf('mean share                           = %7.3f %7.3f \n',nanmean(Share)             ,all_data_moments(5)); 
fprintf('s.d. share                           = %7.3f %7.3f \n',nanstd(Share)              ,all_data_moments(6));
fprintf('median share                         = %7.3f %7.3f \n',nanmedian(Share)           ,all_data_moments(7)); 
fprintf('75th share                           = %7.3f %7.3f \n',prctile(Share,75)          ,all_data_moments(8)); 
fprintf('95th share                           = %7.3f %7.3f \n',prctile(Share,95)          ,all_data_moments(9));
fprintf('99th share                           = %7.3f %7.3f \n',prctile(Share,99)          ,all_data_moments(10));
fprintf('\n');
fprintf('frac. sales top 0.01 sectors         = %7.3f %7.3f \n',fRtopR1_s                  ,all_data_moments(11));
fprintf('frac. sales top 0.05 sectors         = %7.3f %7.3f \n',fRtopR5_s                  ,all_data_moments(12));
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',fRtopL1_s                  ,all_data_moments(13));
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',fRtopL5_s                  ,all_data_moments(14));
fprintf('\n');
fprintf('number producers p10                 = %7.0f %7.0f \n',prctile(Ndom,10)           ,all_data_moments(15));
fprintf('number producers p25                 = %7.0f %7.0f \n',prctile(Ndom,25)           ,all_data_moments(16));
fprintf('number producers p50                 = %7.0f %7.0f \n',prctile(Ndom,50)           ,all_data_moments(17));
fprintf('number producers p75                 = %7.0f %7.0f \n',prctile(Ndom,75)           ,all_data_moments(18));
fprintf('number producers p90                 = %7.0f %7.0f \n',prctile(Ndom,90)           ,all_data_moments(19));
fprintf('number producers p95                 = %7.0f %7.0f \n',prctile(Ndom,95)           ,all_data_moments(20));
fprintf('\n');
fprintf('frac. value added top 0.01 establ.   = %7.3f %7.3f \n',fRtopY1                    ,all_data_moments(21));
fprintf('frac. value added top 0.05 establ.   = %7.3f %7.3f \n',fRtopY5                    ,all_data_moments(22));
fprintf('frac. WL of (same) top 0.01          = %7.3f %7.3f \n',fRtopL1                    ,all_data_moments(23));
fprintf('frac. WL of (same) top 0.05          = %7.3f %7.3f \n',fRtopL5                    ,all_data_moments(24));
fprintf('\n');
fprintf('agg. fraction exporters              = %7.3f %7.3f \n',agg_fexporters             ,all_data_moments(25));
fprintf('agg. import share                    = %7.3f %7.3f \n',agg_impshare               ,all_data_moments(26));
fprintf('\n');
fprintf('coeff. inverse markup wrt share      = %7.3f %7.3f \n',beta_invmarkup_marketshare ,all_data_moments(27)); 
fprintf('\n');
fprintf('index import share dispers.          = %7.3f %7.3f \n',aweight                    ,all_data_moments(28));
fprintf('Grubel-Lloyd (GL) index              = %7.3f %7.3f \n',intraindex                 ,all_data_moments(29));
fprintf('\n');
fprintf('corr(share exports , share sales)    = %7.3f %7.3f \n',corr_share_exports_share_sales  ,all_data_moments(30));  
fprintf('corr(share imports , share sales)    = %7.3f %7.3f \n',corr_share_imports_share_sales  ,all_data_moments(31));  
fprintf('\n');
fprintf('coeff. share imports wrt sh. sales   = %7.3f %7.3f \n',beta_share_imports_share_sales  ,all_data_moments(32));   
fprintf('coeff. share imports wrt sh, exports = %7.3f %7.3f \n',beta_share_imports_share_exports,all_data_moments(33));   
fprintf('\n');
fprintf('mean   p90-p10                       = %7.3f %7.3f \n',nanmean(Share90-Share10)   ,all_data_moments(34)); 
fprintf('median p90-p10                       = %7.3f %7.3f \n',nanmedian(Share90-Share10) ,all_data_moments(35)); 
fprintf('mean   p75-p25                       = %7.3f %7.3f \n',nanmean(Share75-Share25)   ,all_data_moments(36)); 
fprintf('median p75-p25                       = %7.3f %7.3f \n',nanmedian(Share75-Share25) ,all_data_moments(37)); 
fprintf('mean   sd share                      = %7.3f %7.3f \n',nanmean(Sdshare)           ,all_data_moments(38)); 
fprintf('median sd share                      = %7.3f %7.3f \n',nanmedian(Sdshare)         ,all_data_moments(39));
fprintf('\n');
fprintf('HH inverse p10                       = %7.3f %7.3f \n',prctile(IndustryIHH,10)    ,all_data_moments(40));
fprintf('HH inverse p25                       = %7.3f %7.3f \n',prctile(IndustryIHH,25)    ,all_data_moments(41));
fprintf('HH inverse p50                       = %7.3f %7.3f \n',prctile(IndustryIHH,50)    ,all_data_moments(42));
fprintf('HH inverse p75                       = %7.3f %7.3f \n',prctile(IndustryIHH,75)    ,all_data_moments(43));
fprintf('HH inverse p90                       = %7.3f %7.3f \n',prctile(IndustryIHH,90)    ,all_data_moments(44));
fprintf('HH inverse p95                       = %7.3f %7.3f \n',prctile(IndustryIHH,95)    ,all_data_moments(45));
fprintf('\n');
fprintf('top share  p10                       = %7.3f %7.3f \n',prctile(IndustryTop,10)    ,all_data_moments(46)); 
fprintf('top share  p25                       = %7.3f %7.3f \n',prctile(IndustryTop,25)    ,all_data_moments(47));
fprintf('top share  p50                       = %7.3f %7.3f \n',prctile(IndustryTop,50)    ,all_data_moments(48));
fprintf('top share  p75                       = %7.3f %7.3f \n',prctile(IndustryTop,75)    ,all_data_moments(49));
fprintf('top share  p90                       = %7.3f %7.3f \n',prctile(IndustryTop,90)    ,all_data_moments(50));
fprintf('top share  p95                       = %7.3f %7.3f \n',prctile(IndustryTop,95)    ,all_data_moments(51));
fprintf('\n');
fprintf('Trade Elasticity                     = %7.3f %7.3f \n',sigma    , 4);
fprintf('\n');
fprintf('\n');


if 0 %% loss function

loss = (all_data_moments-all_model_moments)./all_data_moments; 
        
% ad hoc weighting
   
loss(25:26) = loss(25:26)*10;   % ensure match export/import data better 
loss(27)    = 0;                % ignore coefficient on sectoral share share
loss(28:29) = loss(28:29)*5;
loss(30:31) = 0;                % ignore regression coefficients
loss(21:22) = loss(21:22)*2;    % worry more about value added concentration
loss(11:12) = loss(11:12)*2; 
loss = [loss; (sigma/4 - 1)]; 
loss(end)   = loss(end)*20;     % lots of weight on trade elasticity


fprintf('\n');
fprintf('\n');
fprintf('mean absolute percentage loss        = %7.3f   \n'    ,nanmean(abs(100*loss)));
fprintf('\n');   
fprintf('\n'); 

end

if 0 %% some other statistics we refer to

beta_imp_HH    
    
disp('aggregate labor share')
agg_lshare=[sum(labor_s)/sum(sales_s)]

disp('average labor share')
avg_lshare = mean(labor./sales)

fprintf('\n')
display('                                         Model    Data')
fprintf('\n')
display('Sectoral markup moments')
fprintf('\n')
fprintf('sectoral markup, p50                   %7.2f %7.2f \n', prctile(muj,50),1.61) % data moments copied from Virgiliu's Jan 2014 email
fprintf('sectoral markup, p75 to p50            %7.2f %7.2f \n', prctile(muj,75)/prctile(muj,50),1.14)
fprintf('sectoral markup, p90 to p50            %7.2f %7.2f \n', prctile(muj,90)/prctile(muj,50),1.43)
fprintf('sectoral markup, p95 to p50            %7.2f %7.2f \n', prctile(muj,95)/prctile(muj,50),1.82)
fprintf('sectoral markup, p99 to p50            %7.2f %7.2f \n', prctile(muj,99)/prctile(muj,50),3.16)
fprintf('\n')
    
 
end
