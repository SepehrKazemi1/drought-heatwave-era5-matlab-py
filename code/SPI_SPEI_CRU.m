%% SPI & (P-PET) based SPEI — 4 separate region figures (NW/NE/SW/SE)
clear; clc;

fPre = "C:\Users\Asus\Downloads\cru_ts4.08.1901.2023.pre.dat.nc";
fPet = "C:\Users\Asus\Downloads\cru_ts4.08.1901.2023.pet.dat.nc";

refStart   = datetime(1980,1,1);
refEnd     = datetime(2010,12,31);
studyStart = datetime(1980,1,1);
studyEnd   = datetime(2020,12,31);
scales     = [3 6 12];
minRefPerMonth = 5;
set(groot,'defaultAxesFontWeight','bold','defaultTextFontWeight','bold','defaultLegendFontWeight','bold','defaultAxesLineWidth',1.2,'defaultAxesFontSize',10);
set(groot,'defaultLineLineWidth',1.2);

% region names (colors optional)
REGNAM = {'NW','NE','SW','SE'};
REGCOL = [0.4940 0.1840 0.5560;
          0.0000 0.4470 0.7410;
          0.8500 0.3250 0.0980;
          0.4660 0.6740 0.1880];

pre  = ncread(fPre, 'pre');       % [lon x lat x time], mm/month
lat  = ncread(fPre, 'lat');
lon  = ncread(fPre, 'lon');
time = ncread(fPre, 'time');
tun  = ncreadatt(fPre, 'time', 'units');

% time decode
basestr = strtrim(extractAfter(tun,'since'));
try, tref = datetime(basestr,'InputFormat','yyyy-MM-dd HH:mm:ss');
catch, tref = datetime(basestr,'InputFormat','yyyy-MM-dd'); end
if contains(lower(tun),'day')
    dates = tref + days(double(time));
elseif contains(lower(tun),'month')
    dates = dateshift(tref,'start','month') + calmonths(double(time));
else
    error('Unexpected time units: %s', tun);
end
pre(pre<-1e10) = NaN;

pet = ncread(fPet,'pet'); pet(pet<-1e10) = NaN;

% study date
tmask = (dates >= studyStart) & (dates <= studyEnd);
dates = dates(tmask);
pre   = pre(:,:,tmask);
pet   = pet(:,:,tmask);

% regions
NW = [32.5 33.5 50.5 52.0];
NE = [32.5 33.5 52.0 53.5];
SW = [31.5 32.5 50.5 52.0];
SE = [31.5 32.5 52.0 53.5];
BOXES = {NW,'NW'; NE,'NE'; SW,'SW'; SE,'SE'};

%  region series
[nlon,nlat,T] = size(pre);
P2 = reshape(pre, nlon*nlat, T);
E2 = reshape(pet, nlon*nlat, T);

regP   = cell(1,4);     % precip
regBal = cell(1,4);     % P - PET

for r = 1:4
    b = BOXES{r,1};
    idx = boxLinearIndex(lon,lat,b);
    regP{r}   = mean(P2(idx,:),1,'omitnan')';
    regBal{r} = (mean(P2(idx,:),1,'omitnan') - mean(E2(idx,:),1,'omitnan'))';
end

%  COMPUTE SPI & SPEI (P-PET via compute_spi_gamma) 
SPI  = cell(4,numel(scales));
SPEI = cell(4,numel(scales));
for r = 1:4
    for s = 1:numel(scales)
        k = scales(s);
        SPI{r,s}  = compute_spi_gamma(regP{r},   dates, k, refStart, refEnd, minRefPerMonth);
        SPEI{r,s} = compute_spi_gamma(regBal{r}, dates, k, refStart, refEnd, minRefPerMonth);
    end
end

%make 4 figs
ylims = [-3 3];

for r = 1:4
    figure('Color','w','Name',sprintf('Isfahan %s — SPI & SPEI (3/6/12)', REGNAM{r}));
    tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    ttxt = title(tl, sprintf('%s — SPI & SPEI (3/6/12)  Baseline %s–%s', ...
        REGNAM{r}, datestr(refStart,'yyyy'), datestr(refEnd,'yyyy')));
    ttxt.FontWeight = 'bold';
    xl = xlabel(tl,'Year'); xl.FontWeight = 'bold';

    ax = gobjects(6,1);
    for s = 1:numel(scales)
        % left column: SPI
        ax(2*s-1) = nexttile; hold on
        plot(dates, SPI{r,s}, 'Color', [0 0 0], 'LineWidth', 1.2);
        yline(0,'k-'); grid on; box on; ylim(ylims);
        ttl = title(sprintf('SPI-%d',scales(s))); ttl.FontWeight = 'bold';
        if s<3, set(gca,'XTickLabel',[]); end
        yl = ylabel('Index'); yl.FontWeight = 'bold';
        set(gca,'FontWeight','bold','LineWidth',1.2);

        % right column: SPEI
        ax(2*s) = nexttile; hold on
        plot(dates, SPEI{r,s}, 'Color', REGCOL(r,:), 'LineWidth', 1.2);
        yline(0,'k-'); grid on; box on; ylim(ylims);
        ttl = title(sprintf('SPEI-%d',scales(s))); ttl.FontWeight = 'bold';
        if s<3, set(gca,'XTickLabel',[]); end
        yl = ylabel('Index'); yl.FontWeight = 'bold';
        set(gca,'FontWeight','bold','LineWidth',1.2);
    end
    linkaxes(ax,'x'); xlim([dates(1) dates(end)]);
end

disp('Done: created 4 separate figures (NW, NE, SW, SE), each with SPI/SPEI at 3/6/12 months, all bold, shared x-axis.')

%%  helper for indices of grid cells be inside a lat/lon box
function idx = boxLinearIndex(lon,lat,box)
% box = [lat_min lat_max lon_min lon_max]
[Lon,Lat] = ndgrid(lon,lat);
M = (Lon>=box(3) & Lon<=box(4) & Lat>=box(1) & Lat<=box(2));
idx = find(M);
end
