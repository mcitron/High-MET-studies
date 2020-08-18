# H/E studies (forked from original High-MET-studies script from Esh)

Macro for producing (many) plots in high-MET ZeroBias samples. The code can be run with

```bash
./compile.sh
./doHOvEStudy.exe <inFile> <outFile> <min pu> <rate above L1 HT> <gen decay in HCAL> <subsystem>
```
once you've done `cmsenv` in a CMSSW release (CMSSW_11_0_2). The options in more details are

* inFile - input file
* outFile - output file
* min pu - only consider events with PU above this value
* rate above L1 HT - bool (0/1) find rate/efficiency on top of L1 HT seed
* gen decay in HCAL - bool (0/1) for signal require decay in HCAL
* subsystem - string can be barrel/endcap/both for which HCAL subsystem you want to look at (see large diff between barrel and endcap)
