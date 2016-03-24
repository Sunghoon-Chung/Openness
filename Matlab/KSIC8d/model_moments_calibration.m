fprintf('...calculating model moments... \n');

%%%%% A. Domestic market share moments   

omegad = omegaH;
omegad = omegaH./repmat(sum(omegaH),N,1); % the market shares of H firms within domestic sectors (NxS)
omegad = omegad(:,sum(omegaH)>0);         % keep only sectors with (strictly) positive share (NxSS as below)

[NN,SS] = size(omegad);                   % NN=N, SS<=S

omegad(omegad==0) = NaN;                  % mimic data structure (use of NaN instead of zeros)

sharedom = omegad;                        % domestic market shares (NxSS)

hh            = nansum(sharedom.^2).*10000;           % sector-level Herfindahl Index (HHI) (1xSS)
highest_share = nanmax(sharedom);                     % highest share in sectors (1xSS)

share         = sharedom(:); share = share(share>0);  % vectorize and keep only positive shares


for s = 1:SS,
    
    share90(s)          = prctile(sharedom(:,s),90);
    share75(s)          = prctile(sharedom(:,s),75);
    share50(s)          = prctile(sharedom(:,s),50);
    share25(s)          = prctile(sharedom(:,s),25);
    share10(s)          = prctile(sharedom(:,s),10);
        
    sdshare(s)          = nanstd(sharedom(:,s));
        
    inversehh(s)        = 1./nansum(sharedom(:,s).^2);      % inverse of HH
        
    topshare(s)         = nanmax(sharedom(:,s));            % top producer's share
    
    CR3(s)              = nansum(sharedom(1:3,s));          % top3 producers' share
    
    ndom(s)             = sum((sharedom(:,s)>0));           % # of positive shares, i.e., # of active producers
                
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

sales_s = sum(sales)';  % sector-level sales
labor_s = sum(labor)';  % sector-level labor inputs

labor = labor(sales>0);
sales = sales(sales>0);

% size distribution of establishments across all sectors in H (data counterpart is value-added. why?)
data = sortrows([sales, labor],1);  % sort [sales, labor] by the first column (ascending order)
fRtopY1  = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); % top 1% establishment share
                                                   % floor: round to the nearest integer less than or equal to element
fRtopY5  = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); % top 5% establishment share
fRtopY10 = sum(data(end-floor(length(data)*0.1)+1:end,1))/sum(data(:,1)); % top 10% establishment share
fRtopY50 = sum(data(end-floor(length(data)*0.5)+1:end,1))/sum(data(:,1)); % top 50% establishment share


% size distribution of sectors (data counterpart is sales)
data  = sortrows([sales_s, labor_s], 1);
fRtopR1_s  = sum(data(end-floor(length(data)*0.01)+1:end,1))/sum(data(:,1)); % top 1% domestic market share
fRtopR5_s  = sum(data(end-floor(length(data)*0.05)+1:end,1))/sum(data(:,1)); % top 5% domestic market share
fRtopR10_s = sum(data(end-floor(length(data)*0.1)+1:end,1))/sum(data(:,1)); % top 10% domestic market share
fRtopR50_s = sum(data(end-floor(length(data)*0.5)+1:end,1))/sum(data(:,1)); % top 50% domestic market share


%%%%% C. Export/Import moments

fexporters     = sum(omegaH>0 & omegaF>0)./sum(omegaH>0); % share of exporters within sectors
agg_fexporters = sum(omegaH(:)>0 & omegaF(:)>0)/sum(omegaH(:)>0); % aggregate exporter share (share of exporters in H)

imports_s      = sum(pF.*yF);
import_shares  = imports_s./sum(pF.*yF+pH.*yH);  % share of imports within sectors
agg_impshare   = sum(pF(:).*yF(:))/sum(pH(:).*yH(:) + pF(:).*yF(:)); % aggregate import share

impshare   = sum(omegaF)';   % market shares of F firms in domestic (H) sectors (1X2000)
expshare   = sum(omegaHs)';  % market shares of H firms in foreign (F) sectors (1X2000)
ttrade     = impshare + expshare;  % total trade shares in sectors of both countries (1X2000)

%%%%% GL index

intraindex = 1 - mean(sj(ttrade>0)'.*abs(impshare(ttrade>0)-expshare(ttrade>0))./(impshare(ttrade>0) + expshare(ttrade>0)))./mean(sj(ttrade>0)');


% Blocks are to make it easier to move them around if required

% Block #1

meanCR3           = nanmean(CR3);
meanHH            = nanmean(HH);
meaninvHH         = nanmean(IndustryIHH);  % simple mean of inverse of HHI(1xSS) ==> a scalar value 
meanmaxshare      = nanmean(Highest_Share);
medianmaxshare    = nanmedian(Highest_Share);          

Block1            = [meanCR3,meanHH,meaninvHH,meanmaxshare,medianmaxshare];

% Block #2

meanshare         = nanmean(Share);
sdshare           = nanstd(Share); 
medianshare       = nanmedian(Share);
ppshare           = prctile(Share,[75,95,99]);         

Block2            = [meanshare,sdshare,medianshare,ppshare];

% Block #3

Block3         = [fRtopR1_s,fRtopR5_s,fRtopR10_s,fRtopR50_s];


% Block #4: percentile number producers 

ppndom            = prctile(Ndom,[10,25,50,75,90]); 

Block4            = [ppndom];

% Block #5

Block5            = [fRtopY1,fRtopY5,fRtopY10,fRtopY50];

% Block #6

Block6            = [agg_fexporters,agg_impshare,intraindex];

% Block #7: CR3

top3share  = prctile(CR3,[10,25,50,75,90]); 

Block7         = [top3share];

% Block #8: percentile HH

pphh  = prctile(HH,[10,25,50,75,90]); 

Block8         = [pphh];

% Block #9: percentile inverse HH

ppindustryihh  = prctile(IndustryIHH,[10,25,50,75,90]); 

Block9         = [ppindustryihh];

% Block #10: percentile top share (across industries)

ppindustrytop  = prctile(IndustryTop,[10,25,50,75,90]); 

Block10         = [ppindustrytop];


all_model_moments = [Block1,Block2,Block3,Block4,Block5,Block6,Block7,Block8,Block9,Block10];


