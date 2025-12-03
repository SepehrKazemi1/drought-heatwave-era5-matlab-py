# Drought–Heatwave Analysis (Central Iran)

>>>> ⚠️ Note: This repo does **not** contain the full project code yet. Some parts are missing for now while the thesis is still in progress.

This repo contains my MSc work on compound drought–heatwave events in the Zayandehrud Basin and Isfahan province (1980–2024).

## What this code does

- Reads ERA5(-Land) data and station time series
- Calculates SPI and SPEI at different time scales
- Derives several heatwave indices (e.g. duration, frequency, intensity)
- Makes basic plots and maps for checking results

The scripts are written in MATLAB, Python, and Google Earth Engine (GEE).

## Why I wrote this

1. I needed a workflow that I can re-run when I fix data issues or change thresholds. A lot of the code is research style meaning it works but it’s not a polished package.
2. I’m using this repo as part of my portfolio to show what I can do with climate and hydrological data.

## Limitations

- Hard-coded for my study area (central Iran).
- Error handling is minimal.
- Some variable names are in Persian because that’s how I wrote them originally 

Ill clean this up more and release the rest of the code as I move forward with my thesis.
