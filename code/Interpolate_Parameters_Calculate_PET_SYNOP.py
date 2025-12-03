#this code is for interplating Synoptic raw data within the 3 days gap and then PET calc
import numpy as np
import pandas as pd
INFILE  = r"urexcel"
INTERP_OUT = r"outputexcel"
PET_OUT = r"petexcel_output"
MONTH_START = "1980-01"    
MONTH_END   = "2023-12"
#interpolate short gaps (consecutive NaNs) per station
MAX_GAP_DAYS = 3
#FAO56 constants
ALBEDO = 0.23
ANG_A, ANG_B = 0.25, 0.50      #Angstrom
KRS = 0.16                     #Hargreaves Rs coefficient
WIND_HEIGHT_M = 10.0           #ffm measured at 10 m (we convert to 2 m)

# pet calculations
def esat(T):  # kPa
    return 0.6108 * np.exp((17.27*T)/(T+237.3))

def daylight_bits(phi, J):
    dr = 1 + 0.033*np.cos(2*np.pi*J/365.0)
    delta = 0.409*np.sin(2*np.pi*J/365.0 - 1.39)
    ws = np.arccos(-np.tan(phi)*np.tan(delta))
    N = 24.0/np.pi * ws
    return dr, delta, ws, N
def Ra_MJ(phi, J):
    dr, delta, ws, _ = daylight_bits(phi, J)
    return (24*60/np.pi)*0.0820*dr*(ws*np.sin(phi)*np.sin(delta)
           + np.cos(phi)*np.cos(delta)*np.sin(ws))

def Rso_clear(Ra, z):
    return (0.75 + 2e-5*z) * Ra

def Rnl_longwave(Rs, Rso, tmax, tmin, ea):
    sigma = 4.903e-9
    tmaxK, tminK = tmax+273.16, tmin+273.16
    Rso_safe = np.where(Rso<=0, np.nan, Rso)
    Rs_Rso = np.clip(Rs/Rso_safe, 0.0, 1.0)
    return sigma*((tmaxK**4 + tminK**4)/2.0) * (0.34 - 0.14*np.sqrt(np.maximum(ea,0))) * (1.35*Rs_Rso - 0.35)

def pressure_kPa(z): return 101.3*((293.0-0.0065*z)/293.0)**5.26
def gamma_kPaC(P):  return 0.000665*P
def delta_kPaC(T):  return 4098.0*(0.6108*np.exp(17.27*T/(T+237.3)))/(T+237.3)**2
def u10_to_u2(u10): return u10 * 4.87 / np.log(67.8*10 - 5.42)

def eto_pm_daily(tmax,tmin,tmean,u2,es,ea,Rn,gamma,delta):
    T=tmean
    num = 0.408*delta*Rn + gamma*(900.0/(T+273.0))*u2*(es-ea)
    den = delta + gamma*(1 + 0.34*u2)
    return num / np.where(den<=0, np.nan, den)

def eto_hargreaves(tmax,tmin,tmean,Ra):
    dT = np.maximum(tmax - tmin, 0)
    return 0.0023*(tmean+17.8)*np.sqrt(dT)*Ra

df = pd.read_excel(INFILE)
df.columns = [str(c).strip() for c in df.columns]

date_col = "data" if "data" in df.columns else ("date" if "date" in df.columns else None)
if not date_col:
    raise SystemExit("No 'data' or 'date' column in file.")
dt = pd.to_datetime(df[date_col], errors="coerce")
if dt.notna().sum() < 0.6*len(dt):
    dt = pd.to_datetime(df[date_col], errors="coerce", dayfirst=True)
df["date"] = dt
df = df[df["date"].notna()].copy()

# numeric fields we wanna interpolate parameters only not pet
numcols = [
    "lat","lon","station_elevation",
    "tmax","tmin","tm",
    "ffm","umax","umin","um","td_m",
    "radglo24","raddir_24","radsw_24","ss24"
]
for c in numcols:
    if c not in df: df[c] = np.nan
df[numcols] = df[numcols].apply(pd.to_numeric, errors="coerce")

#short gaps epr stations
interp_cols = [
    "tmax","tmin","tm","ffm","umax","umin","um","td_m",
    "radglo24","raddir_24","radsw_24","ss24"
]
def interp_group(g):
    g = g.sort_values("date").set_index("date")
    # keep lon and lat and time
    g[interp_cols] = g[interp_cols].interpolate(method="time",
                                                limit=MAX_GAP_DAYS,
                                                limit_direction="both")
    return g.reset_index()

dfi = (df.groupby("station_name", group_keys=False)
         .apply(interp_group)
         .sort_values(["station_name","date"]))

#save interpoled parameters 
dfi.to_excel(INTERP_OUT, index=False)
print("Wrote interpolated parameters", INTERP_OUT)

#Pet calculation according to FAO
dfi["doy"] = dfi["date"].dt.day_of_year
phi = np.deg2rad(dfi["lat"].to_numpy(float))
J   = dfi["doy"].to_numpy(int)
Ra  = Ra_MJ(phi, J)
Rso = Rso_clear(Ra, dfi["station_elevation"].to_numpy(float))
dfi["Ra"], dfi["Rso"] = Ra, Rso
_, _, _, N = daylight_bits(phi, J)
dfi["N"] = pd.Series(N, index=dfi.index, dtype=float)

# Rs  (all in MJ m-2 day-1)
Rs = np.full(len(dfi), np.nan, float)
m1 = dfi["radglo24"].notna()
if m1.any():
    Rs[m1.values] = dfi.loc[m1,"radglo24"].to_numpy(float) * 0.01  # J/cm2/day -> MJ/m2/day
m2 = np.isnan(Rs) & dfi["raddir_24"].notna()
if m2.any():
    Rs[m2] = dfi.loc[m2,"raddir_24"].to_numpy(float) * 1e-6        # J/m2/day -> MJ/m2/day
m3 = np.isnan(Rs) & dfi["radsw_24"].notna()
if m3.any():
    Rs[m3] = dfi.loc[m3,"radsw_24"].to_numpy(float) * 1e-6
m4 = np.isnan(Rs) & dfi["ss24"].notna() & dfi["N"].notna() & dfi["Ra"].notna()
if m4.any():
    ratio = (np.minimum(dfi.loc[m4,"ss24"], dfi.loc[m4,"N"]) / dfi.loc[m4,"N"]).clip(0,1)
    Rs[m4] = (ANG_A + ANG_B * ratio) * dfi.loc[m4,"Ra"]
#HargreavesRs (still PM) as last resort for Rs
tmax_np = dfi["tmax"].to_numpy(float)
tmin_np = dfi["tmin"].to_numpy(float)
dT = np.maximum(tmax_np - tmin_np, 0)
m5 = np.isnan(Rs) & np.isfinite(dT) & np.isfinite(Ra)
if m5.any():
    Rs[m5] = KRS * np.sqrt(dT[m5]) * Ra[m5]

dfi["Rs"]  = Rs
dfi["Rns"] = (1 - ALBEDO) * dfi["Rs"]

# vapour pressures
tmean = dfi["tm"].to_numpy(float)
tmean = np.where(np.isnan(tmean), (tmax_np + tmin_np)/2.0, tmean)
es_tmax = esat(tmax_np)
es_tmin = esat(tmin_np)
es_mean = (es_tmax + es_tmin) / 2.0

ea = np.full(len(dfi), np.nan, float)
rhmax = dfi["umax"].to_numpy(float)
rhmin = dfi["umin"].to_numpy(float)
if np.isfinite(rhmax).any() or np.isfinite(rhmin).any():
    ea = (es_tmin*(rhmax/100.0) + es_tmax*(rhmin/100.0)) / 2.0
if np.isnan(ea).all():
    rhm = dfi["um"].to_numpy(float)
    ea = es_mean * (rhm/100.0)
if np.isnan(ea).all():
    td = dfi["td_m"].to_numpy(float)
    ea = esat(td)

# net radiation n psychrometrics
Rnl = Rnl_longwave(dfi["Rs"].to_numpy(float), dfi["Rso"].to_numpy(float),
                   np.where(np.isnan(tmax_np), tmean, tmax_np),
                   np.where(np.isnan(tmin_np), tmean, tmin_np),
                   ea)
Rn = dfi["Rns"].to_numpy(float) - Rnl
P = pressure_kPa(dfi["station_elevation"].to_numpy(float))
gamma = gamma_kPaC(P)
delta = delta_kPaC(tmean)

# wind 10 m to 2 m
u10 = np.nan_to_num(dfi["ffm"].to_numpy(float), nan=0.0)
u2  = u10_to_u2(u10)

# PM daily ETo
eto_pm = eto_pm_daily(np.where(np.isnan(tmax_np), tmean, tmax_np),
                      np.where(np.isnan(tmin_np), tmean, tmin_np),
                      tmean, u2, es_mean, ea, Rn, gamma, delta)

# Plan B ONLY when PM impossible for that day
eto_hs = eto_hargreaves(tmax_np, tmin_np, tmean, Ra)
need_pm = np.isfinite(Rn) & np.isfinite(ea) & np.isfinite(u2) & np.isfinite(es_mean) & np.isfinite(delta) & np.isfinite(gamma)
ETo = np.where(need_pm, eto_pm, eto_hs)

dfi["ETo_mm_day"] = ETo

# monthly sums with full month headers
dfi["month"] = dfi["date"].dt.to_period("M").astype(str).str.replace("-", "/")
monthly = (dfi.groupby(["station_name","month"], as_index=False)["ETo_mm_day"]
             .sum(min_count=1)
             .rename(columns={"ETo_mm_day": "PET_mm_month"}))

pivot = monthly.pivot(index="station_name", columns="month", values="PET_mm_month")
all_months = pd.period_range(MONTH_START, MONTH_END, freq="M").strftime("%Y/%m")
pivot = pivot.reindex(columns=all_months)   # ensures no missing headings

with pd.ExcelWriter(PET_OUT, engine="openpyxl") as xw:
    dfi.to_excel(xw, index=False, sheet_name="interpolated_parameters")
    dfi[["station_name","date","ETo_mm_day"]].to_excel(xw, index=False, sheet_name="daily_PET")
    pivot.to_excel(xw, sheet_name="monthly_PET_sum")

print("Wrote PET workbook →", PET_OUT)
print("Months in header:", pivot.columns[0], "…", pivot.columns[-1])
