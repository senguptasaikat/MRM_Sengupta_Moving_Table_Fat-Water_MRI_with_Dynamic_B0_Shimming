MRM Sengupta CMT FWMRI with Dynamic B0 shimming
=============================

CMT FWMRI with Dynamic B0 Shimming is a repository for

1. Reconstructing multiecho Golden and Linear angle, radial, continuously moving table whole body MRI scan data,
2. Perform fat/water separations using 2 different algorithms
3. Calculate slice specific whole-body shims for dynamic shimming

Code to download multiecho CMT MRI raw data
from online storage, reconstruct image volumes, and reproduce
figures presented in the journal Magnetic Resonance in Medicine
publication, *Whole-body Continuously Moving Table Fat-Water MRI with Dynamic B0 shimming at 3 Tesla*, is included.  

KNOWN ISSUES
------------
* **2016-04-22** : Data files moved from shutdown cloud storage service to be available as binaries associated with the github.com v1.1 release of this repository; updated download Python script and download CSV index file; removed MATLAB download files that do not work with zipped data files; updated download instructions; removed google analytics beacon on top level README.md file. Users should update to the latest version of the code repository to run the updated data download script.

SETUP
-----
* See [SETUP.md](./SETUP.md) for more details.

LICENSE
-------
* This work is available under the [MIT License](http://opensource.org/licenses/MIT).
* See [LICENSE](./LICENSE) for more details.

CONTRIBUTORS
------------
* See [contributors.txt](./contributors.txt) for more details.

SEE ALSO
--------
* zenodo.org DOI for replication code [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45009.svg)](http://dx.doi.org/10.5281/zenodo.45009)
* data.mendeley.com DOI for replication data [http://dx.doi.org/10.17632/n6s2yx5wsr.1](http://dx.doi.org/10.17632/n6s2yx5wsr.1)
* [https://github.com/welcheb/Reproducible_Research](https://github.com/welcheb/Reproducible_Research)
