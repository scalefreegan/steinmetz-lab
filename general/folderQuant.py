#!/usr/bin/python

import os
import pandas as pd
import numpy as np
import subprocess
import csv
import datetime as dt

# special for remote server with no Xterm
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def calcDateCumSum(df, startDate = "min", endMax = 3650, rangeFreq = "30D", corrF = 1e9):
    end = pd.to_datetime(dt.date.today())
    if startDate == "min":
        start = min(df.loc[:,"date"])
        if (end - start) > pd.Timedelta('%s days' % endMax):
            start = end - pd.Timedelta('%s days' % endMax)
    else:
        start = pd.to_datetime(startDate)

    rng = pd.date_range(start, end, freq = rangeFreq)
    csum = [df.loc[df.loc[:,"date"] <= j,"size"].sum()/corrF for j in rng] # output in corrF units, eg GB
    return pd.DataFrame({"date" : rng, "usage" : csum})

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

    name2email = {}
    with open("/g/steinmetz/brooks/steinmetzUsers.csv") as f:
        next(f)
        reader = csv.reader(f,delimiter=",")
        #print(reader)
        for x,y in reader:
            name2email[x] = y

    proc = subprocess.Popen(["find `pwd` -type f -exec ls -l -k --time-style=long-iso {} \;"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out2 = pd.DataFrame([s.split(None, 7) for s in out.splitlines()],
        columns = ["permissions","nlinks","owner","group","size","date","time","fname"])
    #out2 = out2.convert_objects(convert_numeric=True)
    out2.loc[:,"size"] = pd.to_numeric(out2.loc[:,"size"])
    out2.loc[:,"date"] = pd.to_datetime(out2.loc[:,"date"])

    owners = [i for i in out2["owner"].unique()]
    owners_pass = []
    for i in owners:
        if i not in name2email.keys(): # send email alert to aaron.brooks@embl.de
            msend1 = 'echo "%s" | mailx -v -s "Add /g/steinmetz/ user" aaron.brooks@embl.de'
            m = "You need to add user %s to /g/steinmetz/brooks/steinmetzUsers.csv" % i
            proc = subprocess.Popen([msend1 % m], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
        else:
            if name2email[i] == "NA": # send email alert to aaron.brooks@embl.de
                msend2 = 'echo "%s" | mailx -v -s "Add /g/steinmetz/ user email" aaron.brooks@embl.de'
                m = "You need to add an email address for user %s to /g/steinmetz/brooks/steinmetzUsers.csv" % i
                proc = subprocess.Popen([msend2 % m], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
            else:
                # compile alert for named individual
                owners_pass.append(i)
    owners_pass = ["brooks"] # for testing
    for i in owners_pass:
        try:
            thisdf = out2.loc[out2.loc[:,"owner"]==i,:]

            # total usage GB
            total = round(thisdf.loc[:,"size"].sum() / 1e6, 2)
            median_other = round(out2.groupby("owner").agg({'size' : np.sum}).median()[0] / 1e6, 2)

            # larger than 1 GB
            thisdf_L = thisdf.loc[thisdf.loc[:,"size"]>=1e6,:]
            thisdf_L = thisdf_L.sort_values("size",ascending=False)
            thisL = (thisdf_L["size"]/1e6).map(str) + " " + thisdf_L["fname"]
            lfiles = ("\n").join(j for j in thisL)
            otherdf_L = out2.loc[out2.loc[:,"size"]>=1e6,:]
            median_other_L = otherdf_L.groupby("owner").size().median()

            # cum sum over time plot
            out3 = calcDateCumSum(thisdf, startDate = "min", endMax = 3650, rangeFreq = "1D", corrF = 1e6)
            myplt = out3.plot(x = "date", y = "usage", kind = "area", legend = False, title = "Tier-1 Data Usage: %s" % i)
            myplt.set_xlabel("date")
            myplt.set_ylabel("usage (GB)")
            plt.savefig("/g/steinmetz/brooks/tier1-usage.png", dpi = 150)

            # send mail
            proc = subprocess.Popen(["df -h /g/steinmetz/"], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
            usage = out.splitlines()[2].split()
            statement = 'Dear %s,\n\nHello from the Steinmetz Tier-1 Storage Robowizard!\n\nMy job is to alert you about data usage on the Steinmetz Tier-1 storage drive. Space on this centrally managed service is expensive and limited. Currently our Tier-1 storage drive is %s full! We are using %s of %s of purchased space.\n\nPlease carefully evaluate your usage. Delete unnessary files and archive old projects and files to Tier-2 at /g/tier2/steinmetz/\n\nSave space for others!!\n\nIf you need any assistance, you can send a mail to Aaron Brooks (aaron.brooks@embl.de)\n\nThanks for your help! You will recieve a monthly update about your usage.\n\nYour stats:\n\nA plot with your usage statistics over time is attached\n\nYou are using: %s GBs. The median usage for all other users is %s GBs.\n\nYou have %s files larger than 1 GB. Other users have %s on average.\n\nYour large files (>= 1GB) include: \n\nSize(GB) File\n%s'

            # send mail
            msend3 = 'echo "%s" | mailx -v -s "Tier-1 Usage Report" -a "/g/steinmetz/brooks/tier1-usage.png" %s'
            m = statement % (i, usage[3], usage[1], usage[0], total, median_other, thisdf_L.size, median_other_L, lfiles)
            proc = subprocess.Popen([msend3 % (m,name2email[i])], stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
        except:
            pass
