# High-MET-studies

Macro for producing (many) plots in high-MET ZeroBias samples. The code can be run with

```bash
root -l -b -q doMetStudy.cxx
```

once you've done `cmsenv` in a CMSSW release. This has only been tested with CMSSW\_9\_2\_8, so I recommend using that for the moment.

There's a test option if you want to test additional features/plots. There's a variable `bool isTest` in `doMetStudy`. By default, it's set to `False`, so runs over 1M+ events. If set to `True`, it runs over only 1000.