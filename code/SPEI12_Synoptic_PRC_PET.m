%% This code gives u the SPEI time series from excel files with prc and pet
% excels
clc; clear;close all;

file = "choose_ur_path_1sheetprc_1sheetpet";
Tprc = readtable(file, 'Sheet','prc', 'VariableNamingRule','preserve');
Tpet = readtable(file, 'Sheet','pet', 'VariableNamingRule','preserve');
%time
tprc = normalize_month_time(Tprc{:,1});
tpet = normalize_month_time(Tpet{:,1});
%choose overlap months
[t, iP, iE] = intersect(tprc, tpet);
Tprc = Tprc(iP,:);
Tpet = Tpet(iE,:);

%choose the synop stations
P_esf    = double(Tprc.("Esfahan"));
P_air    = double(Tprc.("Esfahan (Airport)"));
P_shah   = double(Tprc.("Shahrekord"));
PET_esf  = double(Tpet.("Esfahan"));
PET_air  = double(Tpet.("Esfahan (Airport)"));
PET_shah = double(Tpet.("Shahrekord"));
WB_esf  = P_esf  - PET_esf;
WB_air  = P_air  - PET_air;
WB_shah = P_shah - PET_shah;

%12 month sum
WB12_esf  = rolling_sum_12(WB_esf);
WB12_air  = rolling_sum_12(WB_air);
WB12_shah = rolling_sum_12(WB_shah);

%baseline 
baseStart = datetime(1980,1,1); % change to ur desired baseline 
baseEnd   = datetime(2010,12,31);

% function of simple spei
SPEI_esf  = simple_spei12(WB12_esf,  t, baseStart, baseEnd);
SPEI_air  = simple_spei12(WB12_air,  t, baseStart, baseEnd);
SPEI_shah = simple_spei12(WB12_shah, t, baseStart, baseEnd);

% plot of subplots
figure('Color','w','Units','normalized','Position',[0.08 0.1 0.8 0.8])
subplot(3,1,1)
plot(t, SPEI_esf, 'LineWidth', 1)
hold on; yline(0,'k-'); hold off
grid on
title('SPEI-12 — Esfahan (baseline 1980–2010)')
ylim([-3.5 3.5])

subplot(3,1,2)
plot(t, SPEI_air, 'LineWidth', 1)
hold on; yline(0,'k-'); hold off
grid on
title('SPEI-12 — Esfahan (Airport) (baseline 1980–2010)')
ylim([-3.5 3.5])

subplot(3,1,3)
plot(t, SPEI_shah, 'LineWidth', 1)
hold on; yline(0,'k-'); hold off
grid on
title('SPEI-12 — Shahrekord (baseline 1980–2010)')
ylim([-3.5 3.5])
xlabel('Date')

function t = normalize_month_time(v)
    vs = string(v); ok = false;
    fmts = {'yyyy-MM'}; % make sure it meets the synop date format in ur excel
    for f = fmts
        try
            t = datetime(vs,'InputFormat',f{1});
            ok = true;
            break;
        catch
        end
    end
    if ~ok
        t = datetime(vs);
    end
    t = dateshift(t,'start','month');
end

function S = rolling_sum_12(x)
    S = movsum(x, [11 0], 'omitnan');
    cnt = movsum(isfinite(x), [11 0]);
    S(cnt < 12) = NaN;     % to make sure u have 12 months
end

function Z = simple_spei12(D12, t, baseStart, baseEnd)
    Z = nan(size(D12));
    m = month(t);
    for mm = 1:12
        base_idx = (m == mm) & isfinite(D12) & t >= baseStart & t <= baseEnd;
        if nnz(base_idx) < 10
            continue;
        end
        mu = mean(D12(base_idx));
        sd = std(D12(base_idx));
        if sd == 0, continue; end

        all_idx = (m == mm) & isfinite(D12);
        Z(all_idx) = (D12(all_idx) - mu) / sd;
    end
end
