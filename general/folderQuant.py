#!/usr/bin/python

import os
import pandas as pd
import numpy as np
import subprocess


def du(path):
  """disk usage in human readable format (e.g. '2,1GB')"""
  import subprocess
  return float(subprocess.check_output(['du','-s', path]).split()[0].decode('utf-8'))

def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)

def subdirQuant(path):
    import os
    import pandas as pd
    ftype2compress = ["fa","txt","sam"]
    d = {}
    for root, dirs, files in os.walk(path, topdown=True):
        for name in files:
            thispath = os.path.join(root, name)
            thissize = du(thispath)
            thistype = name.split(".")[1] in ftype2compress
            d[thispath] = [thispath,thissize]
    df =  pd.DataFrame.from_dict(d, orient = "index")
    df.columns = ["size","compressable"]
    df = df.sort("size",ascending = False)
    return df

if __name__ == "__main__":

    proc = subprocess.Popen(["find `pwd` -type f -exec ls -l -k --time-style=long-iso {} \; 2> /dev/null"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out2 = pd.DataFrame([s.split(None, ) for s in out.splitlines()],
        columns = ["permissions","nlinks","owner","group","size","date","time","fname"],
        dtype = "str")
    out2["size"] = out2["size"].astype(float)
