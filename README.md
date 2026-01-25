SAPE: Sub-band Aware Power-Equalization for Cellular Networks based on srsRAN 4G
======

[![Build Status](https://github.com/srsran/srsRAN_4G/actions/workflows/ccpp.yml/badge.svg?branch=master)](https://github.com/srsran/srsRAN_4G/actions/workflows/ccpp.yml)
[![CodeQL](https://github.com/srsran/srsRAN_4G/actions/workflows/codeql.yml/badge.svg?branch=master)](https://github.com/srsran/srsRAN_4G/actions/workflows/codeql.yml)
[![Coverity](https://scan.coverity.com/projects/28268/badge.svg)](https://scan.coverity.com/projects/srsran_4g_agpl)

srsRAN is an open source 4G software radio suite developed by [SRS](http://www.srs.io).
See the [srsRAN 4G project pages](https://www.srsran.com) for information, guides and project news.

The srsRAN suite includes:
  * srsUE - a full-stack SDR 4G UE application with prototype 5G features
  * srsENB - a full-stack SDR 4G eNodeB application
  * srsEPC - a light-weight 4G core network implementation with MME, HSS and S/P-GW

For application features, build instructions and user guides see the [srsRAN 4G documentation](https://docs.srsran.com/projects/4g/).

For the documentation and implementation of SAPE on srsRAN, please refer to the [SAPE Documents](./SAPE%20Documents) folder.

For the MATLAB simulation models and results, please refer to the [MATLAB Models](./MATLAB%20Models) folder. These scripts require **MATLAB R2024b** due to updates in dependent MathWorks examples.


The MATLAB Channel dataset (.mat file) can be downloaded from [Google Drive](https://drive.google.com/drive/folders/1b9MiAwl0ipeh0tITQEySzdk-6IzgRc1k?usp=sharing).

For the figures and experimental results used in the paper, please refer to the [Paper Materials](./Paper%20Materials) folder.

## Citation

If you use this codebase, please cite our paper:

Dingyu Zhang, Ish Kumar Jain, and Renjie Zhao. 2026. SAPE: Demystifying Sub-band Aware Power-Equalization for Cellular Networks. In *The 27th International Workshop on Mobile Computing Systems and Applications (HotMobile ’26), February 25–26, 2026, Atlanta, GA, USA*. ACM, New York, NY, USA, 6 pages. https://doi.org/10.1145/3789514.3792036

```bibtex
@inproceedings{zhang2026sape,
  author    = {Dingyu Zhang and Ish Kumar Jain and Renjie Zhao},
  title     = {{SAPE: Demystifying Sub-band Aware Power-Equalization for Cellular Networks}},
  booktitle = {The 27th International Workshop on Mobile Computing Systems and Applications (HotMobile '26)},
  pages     = {1--6},
  year      = {2026},
  address   = {Atlanta, GA, USA},
  publisher = {ACM},
  doi       = {10.1145/3789514.3792036},
  url       = {https://doi.org/10.1145/3789514.3792036}
}
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE_SAPE) for details.

