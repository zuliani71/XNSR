# Calibration datasets

The `.mat` files in this directory are the reference outputs for the nine
datasets distributed in `DataIn`.

- Files in `DataCalib/Triang` use triangular smoothing
  (`SMOOTHING_WIN_TYPE = 'T'`).
- Files in `DataCalib/KonnoOhmachi` use Konno–Ohmachi smoothing
  (`SMOOTHING_WIN_TYPE = 'K'`, `b = 40`).

Both families use the same input signals, output frequency grid and
azimuth/dip grid. All numeric scientific fields are finite.

Running any `Software/test_XN_Cruncher_<dataset>.m` script recreates the
triangular result in `DataOut` and automatically checks its scientific
fields against the matching file in `DataCalib/Triang`. The `_CFGFILE`
scripts demonstrate how to supply an explicit XNSR configuration.

An output can also be compared explicitly after adding `Software` to the
MATLAB path:

```matlab
compare_XNSRCalibration( ...
    fullfile('DataOut','CA04.mat'), ...
    fullfile('DataCalib','Triang','CA04.mat'))
```

Official references can be regenerated after adding `Software` to the
MATLAB path:

```matlab
update_XNSRCalibration("CA04","Triang")
update_XNSRCalibration("CA04","KonnoOhmachi")
```

The updater calls `XN_Cruncher`, writes directly to the appropriate
`DataCalib` subdirectory, and validates the method, parameters, dimensions
and numeric fields. It overwrites the selected official reference, so it is
intended only for deliberate calibration maintenance. Normal tests write to
`DataOut` and use `compare_XNSRCalibration`; they never modify `DataCalib`.

`generate_KonnoOhmachiCalibration` remains available only as a compatibility
wrapper around `update_XNSRCalibration(NAME,"KonnoOhmachi")`.
