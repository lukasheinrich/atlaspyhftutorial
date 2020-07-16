import re
import json
import glob
import subprocess
import shlex

def make_harvest_from_result(result, masses):
    return {
        "CLs": result["CLs_obs"],
        "CLsexp": result["CLs_exp"][2],
        "clsd1s": result["CLs_exp"][1],
        "clsd2s": result["CLs_exp"][0],
        "clsu1s": result["CLs_exp"][3],
        "clsu2s": result["CLs_exp"][4],
        "covqual": 3,
        "dodgycov": 0,
        "excludedXsec": -999007,
        "expectedUpperLimit": -1,
        "expectedUpperLimitMinus1Sig": -1,
        "expectedUpperLimitMinus2Sig": -1,
        "expectedUpperLimitPlus1Sig": -1,
        "expectedUpperLimitPlus2Sig": -1,
        "fID": -1,
        "failedcov": 0,
        "failedfit": 0,
        "failedp0": 0,
        "failedstatus": 0,
        "fitstatus": 0,
        "mn1": masses[2],
        "mn2": masses[1],
        "mode": -1,
        "msb": masses[0],
        "nexp": -1,
        "nofit": 0,
        "p0": 0,
        "p0d1s": -1,
        "p0d2s": -1,
        "p0exp": -1,
        "p0u1s": -1,
        "p0u2s": -1,
        "p1": 0,
        "seed": 0,
        "sigma0": -1,
        "sigma1": -1,
        "upperLimit": -1,
        "upperLimitEstimatedError": -1,
        "xsec": -999007,
    }
    
def harvest_results(regions):
    pattern = re.compile("sbottom_(\d+)_(\d+)_(\d+)")

    dataList = []
    for region in regions:
        harvest = []
        files = "results/region{region}.result.sbottom_*_*_*.json".format(
            region = region,
        )
        for fname in glob.glob(files):
            result = json.load(open(fname))
            m = pattern.search(fname)
            masses = list(map(int, m.groups()))
            # only use 60 GeV
            if masses[2] != 60:
                continue
            harvest.append(make_harvest_from_result(result, masses))
        dataList.append(
            ('region{}'.format(region),harvest)
        )
    return dataList

if __name__ == '__main__':
    main()