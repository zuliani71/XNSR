# XNSR DataOut directory

Test results reproduce the hierarchy used by `DataCalib`:

- `DataOut/Triang`: triangular smoothing results.
- `DataOut/KonnoOhmachi`: Konno-Ohmachi smoothing results.

Each `test_XN_Cruncher_<dataset>_<method>.m` script writes one `.mat` file
to the appropriate subdirectory and compares it with the corresponding
official reference in `DataCalib`.
