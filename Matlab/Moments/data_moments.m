
clc;
clear all;
close all;

load product_dsales_1983;

fprintf('...calculating data moments... \n');  %% '\n': go to the next line


HH = []; Highest_Share = []; 

Share = []; Share90 = []; Share75 = []; Share50 = []; Share25 = []; Share10 = [];
Sdshare = []; Ndom = []; IndustryIHH = []; IndustryTop = []; Sales = []; Sales_s = [];


data = dsales';                      %% each column in data is a sector
     
[NN,SS] = size(data);                %% NN = max # of producers, SS = # of sectors

data(isinf(data))=NaN;               %% replace all "inf" with "nan"

sales        = data;

sales_s      = nansum(sales);            %% column sum excluding NaN cells --> total sales of each sector
sector_sales = repmat(sales_s,NN,1);     %% replicate sales_s by NN rows and 1 columns

sharedom       = data./sector_sales;

hh            = nansum(sharedom.^2);     %% sector level HH
highest_share = max(sharedom);

share = sharedom(:)';                    %% reshape 'share' into a row vector

for s = 1:SS,

    share90(s)          = prctile(sharedom(:,s),90);    %% 90 percentile share for each sector
    share75(s)          = prctile(sharedom(:,s),75);
    share50(s)          = prctile(sharedom(:,s),50);
    share25(s)          = prctile(sharedom(:,s),25);
    share10(s)          = prctile(sharedom(:,s),10);

    sdshare(s)          = nanstd(sharedom(:,s));

    inversehh(s)        = 1./nansum(sharedom(:,s).^2);  %% sector level inverse HH

    topshare(s)         = max(sharedom(:,s));

    ndom(s)             = sum(~isnan(sharedom(:,s)));

end
  
    
    %%%%% pool all 4 years into one big data set
    
    HH            = [HH,hh];
    Highest_Share = [Highest_Share,highest_share];
    Share         = [Share, share];
    
    Share90       = [Share90, share90];
    Share75       = [Share75, share75];
    Share50       = [Share50, share50];
    Share25       = [Share25, share25];
    Share10       = [Share10, share10];
    
    Sdshare       = [Sdshare, sdshare];
    
    Ndom          = [Ndom,ndom];
    
    IndustryIHH   = [IndustryIHH,inversehh];
    IndustryTop   = [IndustryTop,topshare];
    
    Sales         = [Sales  , sales(sales>0)'];
    Sales_s       = [Sales_s, sales_s];

    clear share90 share75 share50 share25 share10 sdshare ndom inversehh topshare


% Blocks are to make it easier to move them around if required

% Block #1: Table 1. Panel A. upper left column

meaninvHH      = nanmean(1./HH); 
medianinvHH    = nanmedian(1./HH);
meanmaxshare   = nanmean(Highest_Share);
medianmaxshare = nanmedian(Highest_Share);          

Block1         = [meaninvHH,medianinvHH,meanmaxshare,medianmaxshare];

% Block #2: Table 1. Panel A. lower left column 

meanshare      = nanmean(Share);
sdshare        = nanstd(Share); 
medianshare    = nanmedian(Share);
ppshare        = prctile(Share,[75,95,99]);         

Block2         = [meanshare,sdshare,medianshare,ppshare];

% Block #3: Table 1. Panel A. upper right column 

Sales_s        = sortrows(Sales_s',1);      % sort rows based on the first column
    
fRtopR1_s      = sum(Sales_s(end-floor(length(Sales_s)*0.01)+1:end,1))/sum(Sales_s(:,1));
fRtopR5_s      = sum(Sales_s(end-floor(length(Sales_s)*0.05)+1:end,1))/sum(Sales_s(:,1));

fRtopL1_s      = 0;
fRtopL5_s      = 0;                             

Block3         = [fRtopR1_s,fRtopR5_s,fRtopL1_s,fRtopL5_s];

% Block #4: percentile number producers 

ppndom         = prctile(Ndom,[10,25,50,75,90,95]); 

Block4         = [ppndom];

% Block #5
    
fRtopY1        = 0; 
fRtopY5        = 0; 

fRtopL1        = 0; 
fRtopL5        = 0;                             
 
Block5         = [fRtopY1,fRtopY5,fRtopL1,fRtopL5];

% Block #6 

agg_fexporters = 0; 
agg_impshare   = 0;

aweight        = 0;
intraindex     = 0;

corr_export_domestic_shares = 0;
corr_import_domestic_shares = 0; 

beta_import_domest_shares   = 0;
beta_import_export_shares   = 0;

Block6         = [agg_fexporters,agg_impshare,aweight,intraindex,corr_export_domestic_shares,corr_import_domestic_shares,beta_import_domest_shares,beta_import_export_shares];


% Block #7: distribution of subsample

mean90to10     = nanmean(Share90-Share10);
median90to10   = nanmedian(Share90-Share10);
mean75to25     = nanmean(Share75-Share25);
median75to25   = nanmedian(Share75-Share25);
meanSDshare    = nanmean(Sdshare);
medianSDshare  = nanmedian(Sdshare);                

Block7         = [mean90to10,median90to10,mean75to25,median75to25,meanSDshare,medianSDshare];

% Block #8: percentile inverse HH

ppindustryihh  = prctile(IndustryIHH,[10,25,50,75,90,95]); 

Block8         = [ppindustryihh];

% Block #9: percentile top share (across industries)

ppindustrytop  = prctile(IndustryTop,[10,25,50,75,90,95]); 

Block9         = [ppindustrytop];


all_data_moments = [Block1,Block2,Block3,Block4,Block5,Block6,Block7,Block8,Block9]';

%  save data moments for future use (comment out as required)

savefile = 'saved_data_moments.mat';
save(savefile,'all_data_moments');


%%%%% report moments


fprintf('\n');
fprintf('Data Moments \n');
fprintf('\n');
fprintf('mean HH inverse                      = %7.3f \n',nanmean(1./HH));  
fprintf('median HH inverse                    = %7.3f \n',nanmedian(1./HH));
fprintf('mean highest share                   = %7.3f \n',nanmean(Highest_Share));  
fprintf('median highest share                 = %7.3f \n',nanmedian(Highest_Share));
fprintf('\n');
fprintf('mean share                           = %7.3f \n',nanmean(Share)); 
fprintf('s.d. share                           = %7.3f \n',nanstd(Share));
fprintf('median share                         = %7.3f \n',nanmedian(Share)); 
fprintf('75th share                           = %7.3f \n',prctile(Share,75)); 
fprintf('95th share                           = %7.3f \n',prctile(Share,95));
fprintf('99th share                           = %7.3f \n',prctile(Share,99));
fprintf('\n');
fprintf('frac. sales top 0.01 sectors         = %7.3f \n',fRtopR1_s); 
fprintf('frac. sales top 0.05 sectors         = %7.3f \n',fRtopR5_s); 
fprintf('frac. WL of (same) top 0.01          = %7.3f \n',fRtopL1_s); 
fprintf('frac. WL of (same) top 0.05          = %7.3f \n',fRtopL5_s);
fprintf('\n');
fprintf('number producers p10                 = %7.0f \n',prctile(Ndom,10));
fprintf('number producers p25                 = %7.0f \n',prctile(Ndom,25));
fprintf('number producers p50                 = %7.0f \n',prctile(Ndom,50));
fprintf('number producers p75                 = %7.0f \n',prctile(Ndom,75));
fprintf('number producers p90                 = %7.0f \n',prctile(Ndom,90));
fprintf('number producers p95                 = %7.0f \n',prctile(Ndom,95));
fprintf('\n');
fprintf('frac. value added top 0.01 establ.   = %7.3f \n',fRtopY1);
fprintf('frac. value added top 0.05 establ.   = %7.3f \n',fRtopY5);
fprintf('frac. WL of (same) top 0.01          = %7.3f \n',fRtopL1); 
fprintf('frac. WL of (same) top 0.05          = %7.3f \n',fRtopL5);
fprintf('\n');
fprintf('agg. fraction exporters              = %7.3f \n',agg_fexporters); 
fprintf('agg. import share                    = %7.3f \n',agg_impshare); 
fprintf('\n');
fprintf('index import share dispers.          = %7.3f \n',1-aweight); % we switched definition to 1- what we had before
fprintf('Grubel-Lloyd (GL) index              = %7.3f \n',intraindex); 
fprintf('\n');
fprintf('corr(export share,domestic share)    = %7.3f \n',corr_export_domestic_shares); 
fprintf('corr(import share,domestic share)    = %7.3f \n',corr_import_domestic_shares); 
fprintf('\n');
fprintf('coeff. import share wrt domest share = %7.3f \n',beta_import_domest_shares);
fprintf('coeff. import share wrt export share = %7.3f \n',beta_import_export_shares);
fprintf('\n');
fprintf('mean   p90-p10                       = %7.3f \n',nanmean(Share90-Share10)); 
fprintf('median p90-p10                       = %7.3f \n',nanmedian(Share90-Share10)); 
fprintf('mean   p75-p25                       = %7.3f \n',nanmean(Share75-Share25)); 
fprintf('median p75-p25                       = %7.3f \n',nanmedian(Share75-Share25)); 
fprintf('mean   sd share                      = %7.3f \n',nanmean(Sdshare)); 
fprintf('median sd share                      = %7.3f \n',nanmedian(Sdshare));
fprintf('\n');
fprintf('HH inverse p10                       = %7.3f \n',prctile(IndustryIHH,10));
fprintf('HH inverse p25                       = %7.3f \n',prctile(IndustryIHH,25));
fprintf('HH inverse p50                       = %7.3f \n',prctile(IndustryIHH,50));
fprintf('HH inverse p75                       = %7.3f \n',prctile(IndustryIHH,75));
fprintf('HH inverse p90                       = %7.3f \n',prctile(IndustryIHH,90));
fprintf('HH inverse p95                       = %7.3f \n',prctile(IndustryIHH,95));
fprintf('\n');
fprintf('top share  p10                       = %7.3f \n',prctile(IndustryTop,10)); 
fprintf('top share  p25                       = %7.3f \n',prctile(IndustryTop,25));
fprintf('top share  p50                       = %7.3f \n',prctile(IndustryTop,50));
fprintf('top share  p75                       = %7.3f \n',prctile(IndustryTop,75));
fprintf('top share  p90                       = %7.3f \n',prctile(IndustryTop,90));
fprintf('top share  p95                       = %7.3f \n',prctile(IndustryTop,95));



           