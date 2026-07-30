# XNSR
# Introduction:
**XNSR** is a **Matlab-based** software designed to compute spectral ratios from single seismic stations (typically Horizontal-to-Vertical Spectral Ratio **HVSR**), introducing a novel and more comprehensive approach. The standard HVSR technique assumes that the maximum seismic response occurs along the horizontal plane, while the minimum lies along the vertical axis. However:
- this assumption may fail in the presence of body waves or significant lateral discontinuities;
- XNSR explores the **3D space** to find any spatial plane where the spectral ratio reaches a maximum or minimum (hence the acronym ma**X** to mi**N** **S**pectral **R**atio – **XNSR**).

![image](/Images/XNSR.full.png)

XNSR uses 3D signal rotations to:<br>
- Identify the plane where the maximum/minimum spectral ratio occurs.
- The process is divided into 10 key steps, including:
	1.	Reading the dataset (SAC or TXT)
	2.	Filtering and tapering
	3.	Signal segmentation
	4.	FFT computation
	5.	Azimuth and dip rotation
	6.	Spectrum computation
	7.	Konno-Omachi filtering
	8.	Ratio calculation (XY/Z → XN)
	9.	Averaging and standard deviation
	10.	2D/3D plotting

The implementation in Matlab is notable for its efficiency and adaptability:
- it combines **parallel processing** and **vectorized operations** to significantly reduce computation time (few minutes for full grid searches);
- the code is structured to handle signal rotations in 3D (**azimuth** and **dip**), producing a high-resolution spectral ratio surface;
- the most CPU-intensive steps (e.g., **Konno-Omachi** filtering) are optimized using Matlab’s **Parallel Computing Toolbox**;
- the software includes **interactive plotting tools** for visual exploration of results, allowing the user to click on specific azimuth-dip combinations to view detailed spectra and statistics.

## Structure:
The XNSR distribution is made of the following dirs and files:
- **_Software_**: it includes all the matlab scripts and functions to perform the XNSR analysis. The main scripts are:
    - **XN_Cruncher.m**: this is the main script, run as it is to have some infos about the usage. It accepts three input parameters:
        - 1st: it is a list of three input source files containig the East-West (E-W), North-South (N-S) and Up (U) components recorded by a seisimometer (both SAC and TXT files are allowed);  
        - 2nd: it is the .mat output file where the results of the XNSR analysis is saved;
        - 3rd: it is a config file including the parameters used to tune the behavior  of the script (see also the details about the **_Cfg_** dir).
    The output is a matlab variable including all the calculus details. The variable is saved into a .mat file that can be used for further elaborations.
    - **XN_plotmatdata.m**: it plots results saved by XN_Cruncher.m. The
      default representation is 2D; an optional argument selects the 3D
      representation. A GUI lets you select a saved XNSR `.mat` file.
    - **test_XN_Cruncher_[TYPE].m**: a list of matlab scripts that can be run to automatically to test one of the datasets included in the **_DataIn_** dir. [TYPE] is one of the daset types available (e.g. Polignano, Lorca, etc). Start  **test_XN_Cruncher_Full.m** if you want to run XNSR using all the datasets included inside the **_DataIn_** dir. The ouputs will be deployed automatically inside the **_DataOut_** dir. This test scripts are very useful to understand the usage of **XN_Cruncher.m**. 
- **_Cfg_**: XN_Cruncher.m needs different parameters which can be set inside a config (cfg) file. **_Cfg_** includes some examples of config file. XN_Cruncher.m accept a 3rd input parameter which is the cfg file. If the cfg file is excluded the script will use a set of defaults parameters embedded inside the code.
- **_Datain_**: in this dir you can find groups of three files belonging to different recordings (sites and experiments). Each group is made of the East-West (E-W), North-South (N-S) and Up (U) components recorded by a seisimometer. The dir includes both SAC and TXT files.
- **_DataOut_**: this dir should include just a README.md file when you clone the XNSR repo on your computer. It is used to automatically save all the .mat output files poroduced by the **test_XN_Cruncher_[TYPE].m** scripts
- **_DataCalib_**: this dir includes pre-elaborated `.mat` outputs produced
  by **XN_Cruncher.m** over the datasets inside **_DataIn_**.
  **_DataCalib/Triang_** contains the references generated with triangular
  smoothing, while **_DataCalib/KonnoOhmachi_** contains the corresponding
  references generated with Konno–Ohmachi smoothing.
- **_Images_**: it includes images used inside the README.md file.
- **_README.md_**: this readme file.

## First run:
- clone the GITHUB XNSR repo on your computer;
- run matlab and add to your matlabpath the XNSR **_Software_** dir;
- do not add external copies of `fft2ft` or `ft2fft`: XNSR includes a
  versioned copy of the common `+spectral` library inside **_Software_**.
  Adding only **_Software_** makes both the legacy entry points and the
  namespaced functions (`spectral.*`) available;
- check that the **_DataOut_** dir just inlcudes the README.md file, otherwise delete all the files included in this dir, except README.md;
- run **test_XN_Cruncher_Full.m**; it will take some minutes depending on
  your hardware. Every dataset script writes its result to **_DataOut_** and
  automatically compares the scientific fields with the corresponding
  reference in **_DataCalib/Triang_**. At the end **_DataOut_** should be
  populated with the following `.mat` files:
    - CA04.mat
    - Edificio_Dorando.mat
    - Ferrara.mat
    - Lorca_001.mat
    - Lorca_002.mat
    - Polignano.mat
    - ValMontanaia_001.mat
    - sanGiuliano_001.mat
    - sanGiuliano_002.mat
- run **XN_plotmatdata.m** and use it to compare the results included in
  **_DataOut_** against the appropriate pre-elaborated analysis in
  **_DataCalib/Triang_** or **_DataCalib/KonnoOhmachi_**.
- Enjoy ;-)
<br>

## Spectral utility dependency

XNSR vendors the common MATLAB spectral utilities in:

```text
Software/+spectral
```

The historical functions `fft2ft` and `ft2fft` remain in `Software` as
compatibility entry points. They delegate the numerical operations to the
vendored package, so a fresh clone is reproducible and does not depend on a
user-specific MATLAB path. To verify the active copy after adding `Software`:

```matlab
which fft2ft -all
which ft2fft -all
which spectral.halfSpectrum
```

The XNSR paths should be listed first. Regression tests are available in:

```matlab
run('Software/test_fft2ft.m')
run('Software/test_triangFilter.m')
run('Software/test_hv_konno.m')
run('Software/test_KonnoOhmachiFilter.m')
run('Software/test_XN_Cruncher_CA04_Konno.m')
```

The scripts named `test_XN_Cruncher_<dataset>.m` are reproducible examples
for the nine bundled experiments. They read the three components from
`DataIn`, write `<dataset>.mat` to `DataOut`, and call
`compare_XNSRCalibration` against `DataCalib/Triang`. The Konno–Ohmachi
CA04 test uses `DataCalib/KonnoOhmachi/CA04.mat`. Custom parameter files
can be passed directly as the third argument of `XN_Cruncher`.

An existing output can also be checked manually:

```matlab
compare_XNSRCalibration( ...
    fullfile('DataOut','CA04.mat'), ...
    fullfile('DataCalib','Triang','CA04.mat'))
```

The comparison covers H/V ratios, frequency vectors, standard deviations,
maxima, orientation vectors and input time series. Plot handles and other
session-specific values are intentionally excluded.

`XN_plotmatdata` uses a 2D representation by default. The initial folder
and representation can be selected independently:

```matlab
XN_plotmatdata()
XN_plotmatdata("3D")
XN_plotmatdata(fullfile(pwd,"DataCalib","Triang"))
XN_plotmatdata(fullfile(pwd,"DataCalib","KonnoOhmachi"),"3D")
```

`hv_konno` is a standalone simplified H/V utility. It uses the shared
one-sided FFT normalization and current Konno–Ohmachi filter, but it is not
called by `XN_Cruncher`, whose workflow additionally performs segmentation,
rotations, statistics and parallel page-wise processing.

Official calibration references are maintained with the general updater:

```matlab
update_XNSRCalibration("CA04","Triang")
update_XNSRCalibration("CA04","KonnoOhmachi")
```

`update_XNSRCalibration` always delegates the scientific calculation to
`XN_Cruncher`. It uses `Cfg/XN_Cruncher_Triang.cfg` for triangular
smoothing and `Cfg/XN_Cruncher_Konno.cfg` for Konno–Ohmachi smoothing,
validates the result, and writes it to `DataCalib/Triang` or
`DataCalib/KonnoOhmachi`. This maintenance function intentionally replaces
the official reference and must not be used as a regression test. The old
`generate_KonnoOhmachiCalibration` entry point is retained only as a
compatibility wrapper.

When a percentage-based triangular window is narrower than the FFT bin
spacing and contains no samples, `triangFilter` uses the closest available
frequency bin. This prevents undefined `0/0` outputs without changing any
frequency whose original triangular window already contained samples.

The Konno-Ohmachi implementation precomputes its normalized smoothing
weights once and reuses them inside the page-wise `parfor`. Filtering is
performed as a matrix product, avoiding the large temporary arrays produced
by the historical `repmat`/`accumarray` implementation while preserving its
numerical result.

`test_XN_Cruncher_CA04_Konno.m` runs the complete 600-second CA04 dataset
with the full azimuth/dip grid, Konno-Ohmachi smoothing and page-wise
`parfor`. It writes `DataOut/CA04_Konno_optimized.mat`, leaving the
historical `DataOut/CA04.mat` untouched.

## References:
- http://dx.doi.org/10.13140/RG.2.2.14803.81443<br>
- http://dx.doi.org/10.13140/RG.2.2.36283.23847<br>
