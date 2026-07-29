# Calibration datasets

The `.mat` files in this directory are the reference outputs for the nine
datasets distributed in `DataIn`.

- Files in `DataCalib` use triangular smoothing (`SMOOTHING_WIN_TYPE = 'T'`).
- Files in `DataCalib/KonnoOhmachi` use Konno–Ohmachi smoothing
  (`SMOOTHING_WIN_TYPE = 'K'`, `b = 40`).

Both families use the same input signals, output frequency grid and
azimuth/dip grid. All numeric scientific fields are finite.

To regenerate one Konno–Ohmachi reference after adding `Software` to the
MATLAB path:

```matlab
generate_KonnoOhmachiCalibration("CA04")
```

The generator writes directly to `DataCalib/KonnoOhmachi` and validates the
method, parameters, dimensions and numeric fields before completing.
