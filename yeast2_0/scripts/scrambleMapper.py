#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import Bio
import datetime
from itertools import repeat
from functools import partial
import multiprocessing as mp
import numpy as np
from pyfaidx import Fasta
import pyfaidx
import pandas as pd
import pybedtools
from pymummer import coords_file, alignment, nucmer
import re
import subprocess
import os
import tempfile
from tqdm import tqdm
import urllib
import uuid
import warnings

from importlib.machinery import SourceFileLoader

def readFASTA(x, splitKey = None):
    """
    Is sequence file? Load from file if so. File should be FASTA format
    Use pyfasta
    """

    if type(x) is not str:
        raise TypeError("input must be type str. filename or sequence")
    if os.path.isfile(x):
        tmp_o = Fasta(x, key_function=lambda key: key.split()[0])
        if (splitKey is None):
            o = tmp_o
        else:
            o = { i.split(splitKey)[0] : tmp_o[i] for i in tmp_o.keys() }
    else:
        o = x
    return o

def readFASTA_SeqIO(x):
    """
    Is sequence file? Load from file if so. File should be FASTA format
    Use SeqIO
    """

    o = []
    if type(x) is list:
        for idx, i in enumerate(x):
            if os.path.isfile(i):
                with open(i, "r") as f:
                    d = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
                    o.append(d)
            else:
                o.append(x)
    elif os.path.isfile(str(x)):
        with open(x, "r") as f:
            o = SeqIO.to_dict(SeqIO.parse(f, "fasta"))
    elif isinstance(x,dict):
        #proper offset
        if isinstance([i for i in x.values()][0], SeqIO.SeqRecord):
            o = x
        else:
            raise TypeError("input must be type filename or SeqIO.SeqRecord)")
    elif  isinstance(x, SeqIO.SeqRecord):
        o = x
    else:
        raise TypeError("input must be type filename or SeqIO.SeqRecord)")
    return o

def writeFASTA(seq, file):

    with open(file, "w") as f:
        if isinstance(seq, pyfaidx.FastaRecord) or isinstance(seq, pyfaidx.Sequence):
            f.write(">" + seq.name + "\n")
            f.write(str(seq) + "\n")
        elif isinstance(seq, pyfaidx.Fasta):
            for i in seq.keys():
                f.write(">" + seq[i].name + "\n")
                f.write(str(seq[i]) + "\n")
        else:
            print("Not pyfaidx. Cannot write at this time.")
    return None

def readGFF(gff):

    gff_cnames = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
    ]
    gff_types = {
        'seqname':str,
        'source':str,
        'feature':str,
        'start':np.int32,
        'end':np.int32,
        'score':str,
        'strand':str,
        'frame':str,
        'attribute':str
    }
    with open(gff, mode = 'r') as f:
        data = [line.strip().split("\t") for line in f if len(line.strip().split("\t")) == 9]
        df = pd.DataFrame.from_records(data, columns = gff_cnames)
        df[['start','end']] = df[['start','end']].apply(pd.to_numeric, errors='ignore')
        df.set_index(["seqname","start","end"], inplace=True)
    return df

def writeGFF(gff,outfile):
    # reset index
    gff = gff.reset_index(level=["seqname", "start", "end"])
    gff = gff.reindex(columns=["seqname",
                        "source",
                        "feature",
                        "start",
                        "end",
                        "score",
                        "strand",
                        "frame",
                        "attribute"
                        ])
    with open(outfile,'w') as f:
        gff.to_csv(f, index=False, sep="\t",header=False)

def getSeq(seq,chr,start,end):
    if os.path.isfile(str(seq)):
        seq = readFASTA_SeqIO(seq)
        #proper offset
        if chr in seq.keys():
            return seq[chr][start-1:end-1]
        else:
            warnings.warn(('Could not find sequence {} in supplied reference fasta file!!!'
                           ' Dropping gff entry').format(chr))
            return None
    elif isinstance(seq,dict):
        #proper offset
        if isinstance([i for i in seq.values()][0], Bio.SeqRecord.SeqRecord):
            if chr in seq.keys():
                return seq[chr][start-1:end-1]
            else:
                warnings.warn(('Could not find sequence {} in supplied reference fasta file!!!'
                               ' Dropping gff entry').format(chr))
                return None
    else:
        raise TypeError("input must be Bio.SeqRecord.SeqRecord")

def runBlat(ref, query, blatPath):
    """
    Align sequences with Blat
    """
    if os.path.isfile(str(ref)):
        reference_file = ref
        reference_name = ref
    elif isinstance(ref,Bio.SeqRecord.SeqRecord):
        reference_file = tempfile.NamedTemporaryFile(mode='w')
        reference_name = reference_file.name
        Bio.SeqIO.write(ref, reference_name, "fasta")
    elif isinstance(ref,dict):
        if isinstance([i for i in ref.values()][0], Bio.SeqRecord.SeqRecord):
            reference_file = tempfile.NamedTemporaryFile(mode='w')
            reference_name = reference_file.name
            towrite = [i for i in ref.values()]
            Bio.SeqIO.write(towrite, reference_name, "fasta")
        else:
            raise TypeError("input must be type filename or Bio.SeqRecord.SeqRecord)")
    else:
        raise TypeError("input must be Bio.SeqRecord.SeqRecord")

    if os.path.isfile(str(query)):
        query_file = query
        query_name = query
    elif isinstance(query,Bio.SeqRecord.SeqRecord):
        query_file = tempfile.NamedTemporaryFile(mode='w')
        query_name = query_file.name
        Bio.SeqIO.write(query, query_name, "fasta")
    elif isinstance(query,dict):
        if isinstance([i for i in query.values()][0], Bio.SeqRecord.SeqRecord):
            query_file = tempfile.NamedTemporaryFile(mode='w')
            query_name = query_file.name
            towrite = [i for i in query.values()]
            Bio.SeqIO.write(towrite, query_name, "fasta")
        else:
            raise TypeError("input must be type filename or Bio.SeqRecord.SeqRecord)")
    else:
        raise TypeError("input must be Bio.SeqRecord.SeqRecord")

    results_file = tempfile.NamedTemporaryFile(mode='w', delete=False)
    #print(results_file.name)
    cmds = [
        blatPath,
        reference_name,
        query_name,
        results_file.name
    ]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError("cmds: %s\nstderr:%s\nstdout:%s"
                         % (" ".join(cmds), stderr, stdout))

#     with open(results_file.name,'r') as f:
#         print(f.readlines())

    df = pd.read_csv(results_file.name, sep ="\t", skiprows=[0,1,4], header=None)
    df.columns = (df.iloc[0].map(str) + df.iloc[1].map(str))
    df.columns = [i.replace(" ","").replace("nan","") for i in df.columns]
    df = df.drop([0,1])
    df = df.reset_index(drop=True)
    return(df)

def runMummer(ref, query):
    """
    Align sequences with MUMmer
    """
    useTempRef = False
    useTempQuery = False

    if os.path.isfile(str(ref)):
        reference_file = ref
        reference_name = ref
    elif isinstance(ref,Bio.SeqRecord.SeqRecord):
        useTempRef = True
        reference_file = tempfile.TemporaryFile(mode='w')
        #NamedTemporaryFile(mode='w+b', delete=False)
        # reference_file = os.path.abspath(str(datetime.datetime.now().time())
        #                     + "_ref.fasta")
        reference_file = str(uuid.uuid4())
        Bio.SeqIO.write(ref, reference_file, "fasta")
    elif isinstance(ref,dict):
        if isinstance([i for i in ref.values()][0], Bio.SeqRecord.SeqRecord):
            useTempRef = True
            #reference_file = tempfile.TemporaryFile(mode='w+t')
            # reference_file = os.path.abspath(str(datetime.datetime.now().time())
            #                     + "_ref.fasta")
            reference_file = str(uuid.uuid4())
            towrite = [i for i in ref.values()]
            Bio.SeqIO.write(towrite, ref_file, "fasta")
        else:
            raise TypeError("input must be type filename or Bio.SeqRecord.SeqRecord)")
    else:
        raise TypeError("input must be Bio.SeqRecord.SeqRecord")

    if os.path.isfile(str(query)):
        query_file = query
        query_name = query
    elif isinstance(query,Bio.SeqRecord.SeqRecord):
        useTempQuery = True
        #query_file = tempfile.TemporaryFile(mode='w+t')
        # query_file = os.path.abspath(str(datetime.datetime.now().time())
        #                     + "_query.fasta")
        query_file = str(uuid.uuid4())
        Bio.SeqIO.write(query, query_file, "fasta")
    elif isinstance(query,dict):
        if isinstance([i for i in query.values()][0], Bio.SeqRecord.SeqRecord):
            useTempQuery = True
            #reference_file = tempfile.TemporaryFile(mode='w+t')
            # query_file = os.path.abspath(str(datetime.datetime.now().time())
            #                     + "_query.fasta")
            query_file = str(uuid.uuid4())
            towrite = [i for i in query.values()]
            Bio.SeqIO.write(towrite, query_file, "fasta")
        else:
            raise TypeError("input must be type filename or Bio.SeqRecord.SeqRecord)")
    else:
        raise TypeError("input must be Bio.SeqRecord.SeqRecord")
    # results_file = os.path.abspath(str(datetime.datetime.now().time())
    #                         + "_results")
    results_file = str(uuid.uuid4())
    runner = nucmer.Runner(reference_file, query_file, results_file)
    runner.run()
    file_reader = coords_file.reader(results_file)
    #alignments = [coord for coord in file_reader if not coord.is_self_hit()] #Remove self hits
    alignments = [coord for coord in file_reader]
    # print(alignments)
    # clean up
    if useTempRef:
        os.remove(reference_file)
    if useTempQuery:
        os.remove(query_file)
    os.remove(results_file)
    # print(reference_file)
    # print(query_file)
    # print(results_file)

    col_labels = ["S1", "E1", "S2", "E2", "LEN_1", "LEN_2", "%_IDY",
        "LEN_R", "LEN_Q", "FRM", "NAME_R", "NAME_Q"]
    if alignments:
        #print([str(i).split("\t") for i in alignments])
        df = pd.DataFrame.from_records([str(i).split("\t") for i in alignments],
                columns=col_labels)
    else:
        df = pd.DataFrame([],columns=col_labels)

    return df

def gffMod(df, hits, doReps=False, repNum=1, use=["blat", "mummer"][0]):

    def add2id(matchobj,i):
        x = str(matchobj.group(0).rstrip(";"))+ "-" + str(i) + ";"
        return x
    out = pd.DataFrame.copy(df)
    if use == "mummer":
        seqname = hits["NAME_Q"]
        start_s = int(hits["S2"])
        end_s = int(hits["E2"])+1
        thisstrand = out.iloc[0]["strand"]
    elif use == "blat":
        seqname = str(hits["Qname"])
        start_s = int(hits["Qstart"])+1
        end_s = int(hits["Qend"])+1
        thisstrand = str(hits["strand"])
    else:
        seqname = hits["seqname"]
        start_s = int(hits["start"])
        end_s = int(hits["end"])
        thisstrand = str(hits["strand"])
    #print(thisstrand)
    oppositestrand = lambda x: "-" if x is "+" else "+"
    # clean up web encoding (affects tRNA genes)
    attribute = urllib.request.unquote(out.iloc[0]["attribute"]).rstrip('"').lstrip('"')
    # do some attributes in yeast gff have two sets of quotes??
    attribute = attribute.rstrip('"').lstrip('"')
    if doReps:
    	attribute = re.sub("ID=.*?;", lambda line: add2id(line,repNum), attribute)
    	attribute = re.sub("Name=.*?;", lambda line: add2id(line,repNum), attribute)
    	attribute = re.sub("gene=.*?;", lambda line: add2id(line,repNum), attribute)
    if end_s < start_s:
        start = end_s
        end = start_s
        strand = oppositestrand(thisstrand)
        #strand = thisstrand
    else:
        start = start_s
        end = end_s
        strand = thisstrand
    out.index = pd.MultiIndex.from_tuples([(seqname,start,end)], names=("seqname","start","end"))
    out.loc[(seqname,start,end),"strand"] = strand
    out.loc[(seqname,start,end),"attribute"] = attribute
    return out

def mapGFF(gff_in, to_fa, from_fa, onlyHighest,
           correct_tRNA, enforceLen, use,
           blatPath, unmapped, useEmergencyMapper = False, full_gff = None):
    def fail(x, unmapped, emergency, full_gff,
             to_fa, from_fa, onlyHighest, correct_tRNA,
             enforceLen, use, blatPath):
        print("Failed to translate gff entry: {}".format(x))
        if emergency:
            print("Trying to map with emergency mapper")
            y = emergencyMapper(x, full_gff, to_fa, from_fa, onlyHighest,
                       correct_tRNA, enforceLen, use,blatPath)
        else:
            y = None
        if y is None:
            if unmapped is not None:
                x.to_csv(unmapped, sep = "\t", mode = "a", header = False)
        return y
    #print(gff_in)
    if type(from_fa) is Bio.SeqRecord.SeqRecord:
        feature_seq = from_fa
    else:
        from_seq = readFASTA_SeqIO(from_fa)
        feature_seq = getSeq(from_seq,gff_in.index.get_level_values('seqname')[0],
                           gff_in.index.get_level_values('start')[0],
                           gff_in.index.get_level_values('end')[0])
    if (use=="blat") and (len(feature_seq) > 5000):
        message = "Feature {} has length {} which is too long for blat. Falling back to use mummer"
        print(message.format(gff_in.attribute[0], len(feature_seq)))
        use = "mummer"
    if (use=="mummer") and (len(feature_seq) <= 1000):
        message = "Feature {} has length {} which is too short for mummer. Falling back to use blat"
        print(message.format(gff_in.attribute[0], len(feature_seq)))
        use = "blat"
    #print(feature_seq)
    if feature_seq is not None:
        to_seq = readFASTA_SeqIO(to_fa)
        if use == "mummer":
            mummerhits = runMummer(feature_seq,to_seq)
            mummerhits["%_IDY"] = pd.to_numeric(mummerhits["%_IDY"])
            #print(mummerhits)
            if mummerhits.shape[0] > 0:
                if mummerhits.apply(lambda x : x["NAME_R"] == x["NAME_Q"], axis=1).any():
                    # exact chr match -- just return that
                    mummerhits = mummerhits[mummerhits.apply(lambda x : x["NAME_R"] == x["NAME_Q"], axis=1)]
                # this is a syn match?
                if onlyHighest:
                    mummerhits = mummerhits[mummerhits["%_IDY"] >= mummerhits.max()["%_IDY"]]
                if correct_tRNA:
                    # correct for tRNA and other annotations that may present mapping problems
                    if gff_in.feature[0]=="tRNA_gene":
                        mummerhits = mummerhits[mummerhits.apply(lambda x : x["NAME_R"] == x["NAME_Q"], axis=1)]
                if enforceLen:
                    # make the ref:query match lengths equal
                    mummerhits = mummerhits[mummerhits.apply(lambda x : x["LEN_1"] == x["LEN_2"], axis=1)]
                #print(mummerhits)
                mummerhits = mummerhits.reset_index()
                if mummerhits.shape[0] > 0:
                    matches = mummerhits
                else:
                    matches = fail(gff_in, unmapped, useEmergencyMapper, full_gff,
                             to_fa, from_fa, onlyHighest, correct_tRNA,
                             enforceLen, use, blatPath)
                    use = None

            else:
                matches = fail(gff_in, unmapped, useEmergencyMapper, full_gff,
                         to_fa, from_fa, onlyHighest, correct_tRNA,
                         enforceLen, use, blatPath)
                use = None
        elif use == "blat":
            blathits = runBlat(feature_seq,to_seq, blatPath)
            blathits["%_IDY"] = pd.to_numeric(blathits["match"])/pd.to_numeric(blathits["Tsize"])
            #print(blathits)
            if blathits.shape[0] > 0:
                if blathits.apply(lambda x : x["Tname"] == x["Qname"], axis=1).any():
                    # exact chr match -- just return that
                    blathits = blathits[blathits.apply(lambda x : x["Tname"] == x["Qname"], axis=1)]
                # this is a syn match?
                if onlyHighest:
                    blathits = blathits[blathits["%_IDY"] >= blathits.max()["%_IDY"]]
                if correct_tRNA:
                    # correct for tRNA and other annotations that may present mapping problems
                    if gff_in.feature[0]=="tRNA_gene":
                        blathits = blathits[blathits.apply(lambda x : x["Tname"] == x["Qname"], axis=1)]
                if enforceLen:
                    # make the ref:query match lengths equal
                    blathits = blathits[blathits.apply(lambda x : x["Tsize"] == x["Qsize"], axis=1)]
                blathits = blathits.reset_index()
                if blathits.shape[0] > 0:
                    matches = blathits
                else:
                    matches = fail(gff_in, unmapped, useEmergencyMapper, full_gff,
                             to_fa, from_fa, onlyHighest, correct_tRNA,
                             enforceLen, use, blatPath)
                    use = None
            else:
                matches = fail(gff_in, unmapped, useEmergencyMapper, full_gff,
                         to_fa, from_fa, onlyHighest, correct_tRNA,
                         enforceLen, use, blatPath)
                use = None
    else:
        matches = fail(gff_in, unmapped, useEmergencyMapper, full_gff,
                 to_fa, from_fa, onlyHighest, correct_tRNA,
                 enforceLen, use, blatPath)
        use = None
    if matches is not None:
        if matches.shape[0] > 1:
            if use == "blat":
                matches["Qend"] = matches["Qend"].astype('int64')
                matches = matches.sort_values("Qend").reset_index()
            elif use == "mummer":
                matches["S2"] = matches["S2"].astype('int64')
                matches = matches.sort_values("S2").reset_index()
            out = gffMod(gff_in, matches.iloc[0], doReps=True, repNum=1, use=use)
            for i in range(1,(matches.shape[0])):
                #print(i)
                newrow = gffMod(gff_in,matches.loc[i], doReps=True, repNum=i+1, use=use)
                out = out.append(newrow)
        else:
            out = gffMod(gff_in, matches.iloc[0], use=use)
    else:
        out = None
    return(out)

def mapMultipleGFF(from_fa, gff_in, to_fa, onlyHighest=True,
           correct_tRNA=True, enforceLen=False, use=["blat","mummer"][0],
           blatPath = "/g/steinmetz/brooks/anaconda/envs/nanopore/bin/blat",
           usemp = True, nproc = 32, unmapped = "unmapped.gff",
           useEmergencyMapper = True):
    with open(unmapped, 'w') as f:
        f.write("# unmapped gff entries\n")
    if gff_in.shape[0] > 0:
        if usemp:
            print("Using multiprocessing")
            pool = mp.Pool(processes=nproc)
            #results = [pool.apply(mapGFF, args=(from_fa,gff_in.iloc[[x]], to_fa)) for x in range(0,gff_in.shape[0])]
            gff_arg = [gff_in.iloc[[x]] for x in range(0,gff_in.shape[0])]
            #results = pool.starmap(mapGFF, zip(gff_arg,repeat(to_fa),repeat(from_fa)))
            fcn = partial(mapGFF, to_fa = to_fa, from_fa = from_fa, onlyHighest = onlyHighest,
                       correct_tRNA = correct_tRNA, enforceLen = enforceLen, use = use,
                       blatPath = blatPath, unmapped = unmapped, full_gff = gff_in,
                       useEmergencyMapper = useEmergencyMapper)
            #print(fcn)
            results = pool.map(fcn, gff_arg)
            #results = pool.map(partial(mapGFF, gff_in=gff_arg), to_fa=to_fa, from_fa=from_fa)
            #results.close()
            #results.join()
            out = pd.concat(results)
            return out
        else:
            print("Not implemented")
            # print("Using tqdm")
            # out = mapGFF(from_fa, gff_in.iloc[[0]], to_fa)
            # if gff_in.shape[0] > 1:
            #     for i in tqdm(range(1,(gff_in.shape[0]))):
            #         #print(i)
            #         newrow = mapGFF(from_fa, gff_in.iloc[[i]], to_fa)
            #         # if newrow is not None:
            #         #     out = out.append(newrow)
            # return out
    else:
        return None

def addSeg(to_fa, seg_file, gff = None, onlyHighest=True,
           correct_tRNA=True, enforceLen=False, use=["blat","mummer"][0],
           blatPath = "/g/steinmetz/brooks/anaconda/envs/nanopore/bin/blat",
           usemp=True, nproc=32, unmapped = "unmapped.gff"):

    if os.path.isfile(seg_file):
        with open(seg_file,'r') as f:
            df = pd.read_csv(f,sep='\t')
        temp_file = tempfile.NamedTemporaryFile(mode='w')
        #temp_file = str(uuid.uuid4())
        with open(temp_file.name,'w') as f2:
            for i in df["number"]:
                f2.write(">"+str(i) + "\n")
                f2.write(str(df.loc[df["number"] == i].seq)+"\n")
        dummy_gff = pd.DataFrame.from_records([{
            "seqname":"dummy",
            "start":1,
            "end":1,
            "source":"ANB",
            "feature":"engineered_region",
            "score":".",
            "strand":"+",
            "frame":".",
            "attribute":"ID=;"
        }],index=("seqname","start","end"))
        thisfrom = readFASTA_SeqIO(temp_file.name)
        #out = pd.DataFrame()
        results = []
        segments = [x for x in thisfrom.keys()]
        if usemp:
            gff_arg = dummy_gff
            for i in range(0,len(segments)):
                ii = segments[i]
                if i == 0:
                    gff_arg.iloc[i].attribute = "ID=" + str(ii) + ";"
                else:
                    gff_arg = gff_arg.append(dummy_gff)
                    gff_arg.iloc[i].attribute = "ID=" + str(ii) + ";"
            gff_arg = [gff_arg.iloc[[x]] for x in range(0,gff_arg.shape[0])]
            from_fa_arg = [thisfrom[s] for s in segments]
            pool = mp.Pool(processes=nproc)
            results = pool.starmap(mapGFF, zip(gff_arg,
                                                repeat(to_fa),
                                                [thisfrom[s] for s in segments],
                                                repeat(onlyHighest),
                                                repeat(correct_tRNA),
                                                repeat(enforceLen),
                                                repeat(use),
                                                repeat(blatPath),
                                                repeat(unmapped),
                                                ))
            # def fcn(i):
            #     result = mapGFF(gff_in = gff_arg[i], to_fa = to_fa, from_fa = from_fa_arg[i],
            #             onlyHighest = onlyHighest,correct_tRNA = correct_tRNA,
            #             enforceLen = enforceLen, use = use,
            #             blatPath = blatPath, unmapped = unmapped)
            #     return(result)
            # results = pool.map(fcn, list(range(0,len(gff_arg))))
            out = pd.concat(results)
        else:
            for i in tqdm(range(0,len(segments))):
                ii = segments[i]
                dummy_gff.iloc[0].attribute = "ID=" + str(ii) + ";"
                thisentry = mapGFF(gff_in = dummy_gff, to_fa = to_fa, from_fa = thisfrom[ii])
                if thisentry is not None:
                    results.append(thisentry)
            out = pd.concat(results)
        if gff is not None:
            out = gff.append(out)
        return out
    else:
        raise IOError("must provide filename for argument segment_table")

def emergencyMapper(x, gff, to_fa, from_fa, onlyHighest,
           correct_tRNA, enforceLen, use,
           blatPath):
    """
    Map entry based on neighboring alignments. Coordinates are relative to
    neighboring alignments. This makes a lot of assumptions
    but is necessary for small fragments that do not have unique
    alignments. Relies on mapping with SCRaMbLE segs.
    Need to do this wrt to segments to other mappings!!
    """

    if "site_specific_recombination_target_region" in list(gff.feature):
        x = x.reset_index()
        seq = x["seqname"][0]
        start = x.start[0]
        end = x.end[0]
        # locate adjacent loxPsym sites
        gff = gff.reset_index()
        loxPsyms = gff.loc[gff.feature == "site_specific_recombination_target_region"]
        loxPsyms = loxPsyms.assign(dstart = abs(loxPsyms["start"] - start),
                                 dend = abs(loxPsyms["end"] - end))
        loxPsyms_start = loxPsyms.sort_values(["dstart"]).reset_index()
        loxPsyms_end = loxPsyms.sort_values(["dend"]).reset_index()

        # start of feature is closest
        if loxPsyms_start["dstart"].iloc[0] <= loxPsyms_end["dend"].iloc[0]:
            if loxPsyms_start["end"].iloc[0] < start:
                leftSym = loxPsyms_start.iloc[[0]]
                # special exception for last loxPsym site
                if (int(leftSym["start"]) >= max(loxPsyms_start["start"])):
                    rightSym = loxPsyms_start.iloc[[0]]
                    rightSym["start"] = max(gff["end"]) + 1
                    rightSym["end"] = max(gff["end"]) + 1
                else:
                    i = 0
                    while (i < loxPsyms_end.shape[0]):
                        while ((loxPsyms_end.iloc[[i]].equals(leftSym)) or
                              (int(loxPsyms_end.iloc[[i]]["start"]) < int(leftSym["start"]))):
                              i = i + 1
                    rightSym = loxPsyms_end.iloc[[i]]
            else:
                rightSym = loxPsyms_start.iloc[[0]]
                # special exception for first loxPsym site
                if (int(rightSym["start"]) <= min(loxPsyms_start["start"])):
                    leftSym = loxPsyms_start.iloc[[0]]
                    leftSym["start"] = 0
                    leftSym["end"] = 0
                else:
                    i = 0
                    while (i < loxPsyms_end.shape[0]):
                        while ((loxPsyms_start.iloc[[i]].equals(rightSym)) or
                              (int(loxPsyms_start.iloc[[i]]["start"]) > int(rightSym["start"]))):
                            i = i + 1
                    leftSym = loxPsyms_start.iloc[[i]]
        # start of feature is closest
        elif loxPsyms_end["dend"].iloc[0] < loxPsyms_start["dstart"].iloc[0]:
            if loxPsyms_end["start"].iloc[0] > end:
                rightSym = loxPsyms_end.iloc[[0]]
                # special exception for first loxPsym site
                if (int(rightSym["start"]) <= min(loxPsyms_start["start"])):
                    leftSym = loxPsyms_start.iloc[[0]]
                    leftSym["start"] = 0
                    leftSym["end"] = 0
                else:
                    i = 0
                    while (i < loxPsyms_end.shape[0]):
                        while ((loxPsyms_start.iloc[[i]].equals(rightSym)) or
                              (int(loxPsyms_start.iloc[[i]]["start"]) > int(rightSym["start"]))):
                            i = i + 1
                    leftSym = loxPsyms_start.iloc[[i]]
            else:
                leftSym = loxPsyms_end.iloc[[0]]
                # special exception for last loxPsym site
                if (int(leftSym["start"]) >= max(loxPsyms_start["start"])):
                    rightSym = loxPsyms_start.iloc[[0]]
                    rightSym["start"] = max(gff["end"]) + 1
                    rightSym["end"] = max(gff["end"]) + 1
                else:
                    i = 0
                    while (i < loxPsyms_end.shape[0]):
                        while ((loxPsyms_end.iloc[[i]].equals(leftSym)) or
                              (int(loxPsyms_end.iloc[[i]]["start"]) < int(leftSym["start"]))):
                            i = i + 1
                    rightSym = loxPsyms_end.iloc[[i]]

        gff_sub = gff.loc[(gff["seqname"]==seq) &
                           (gff["start"] > int(leftSym["end"])) &
                           (gff["end"] < int(rightSym["start"]))]
        gff_sub = gff_sub.assign(dstart = abs(gff_sub["start"] - start),
                                 dend = abs(gff_sub["end"] - end))
        gff_sub = gff_sub.loc[(gff_sub["dstart"]>0) & (gff_sub["dend"]>0)]
        if gff_sub.shape[0] == 0:
            print(x)
            print(gff_sub)
            print("Failed to map features in area")
            toreturn = None
        else:

            if min(gff_sub["dstart"]) < min(gff_sub["dend"]):
                gff_sub = gff_sub.sort_values(["dstart"])
            else:
                gff_sub = gff_sub.sort_values(["dend"])
            i = 0
            closest = None
            while ((closest is None) & (i < gff_sub.shape[0])):
                xi = gff_sub.iloc[[i]].set_index(["seqname","start","end"])
                closest = mapGFF(xi, to_fa, from_fa, onlyHighest,
                           correct_tRNA, enforceLen, use,
                           blatPath, unmapped = None,
                           useEmergencyMapper = False,
                           full_gff = gff)
                i = i + 1
            if closest is None:
                print(x)
                print(gff_sub)
                print("Failed to map features in area")
                toreturn = None
            else:
                thisdistance = int(x["start"]) - int(gff_sub.iloc[(i-1)]["start"])
                closest = closest.reset_index()
                def gffUpdater(this, x = x, thisdistance = thisdistance):
                    this["start"] = int(this["start"]) + thisdistance
                    this["end"] = int(this["end"]) + thisdistance
                    this["feature"] = str(x["feature"][0])
                    this["source"] = str(x["source"][0])
                    this["score"] = str(x["score"][0])
                    this["attribute"] = str(x["attribute"][0])
                    if str(x["strand"][0]) == ".":
                        this["strand"] = str(x["strand"][0])
                    return(this)

                toreturn = closest.apply(gffUpdater, axis=1)
        return(toreturn)
    else:
        print("ERROR: Could not find loxPsym sequences in GFF file. Can't perform emergency mapping")
        return(None)


if __name__ == "__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_fa', help='file path to reference genome, fasta format',
                        required=True, metavar = "REFERENCE FASTA FILE PATH")
    parser.add_argument('--ref_gff', help='file path to reference genome, gff format',
                        required=False, default=None, metavar = "REFERENCE GFF FILE PATH")
    parser.add_argument('--query_fa', help='file path to query genome, fasta format',
                        required=True, metavar = "QUERY FASTA FILE PATH")
    parser.add_argument('--output_gff', help='file path to output gff',
                        required=True, metavar = "OUTPUT GFF FILE PATH")
    parser.add_argument('--annotation_filter', help='list of gff "type" features to exclude',
                        required=False, default=None, nargs="+", metavar = "FEATURE TYPE LIST")
    parser.add_argument('--segment_table', help='segment sequence table file name. columns: number,seq',
                        required=False, default=None, type=str, metavar = "FEATURE TYPE LIST")
    parser.add_argument('--cores', help='number of cores to use',
                        required=False, default=32, type=int)
    parser.add_argument('--use', help='aligner to use, blat or mummer',
                        required=False, default="blat", type=str)
    parser.add_argument('--blatpath', help='path to blat executable',
                        required=False, default="/g/steinmetz/brooks/anaconda/envs/nanopore/bin/blat", type=str)
    parser.add_argument('--unmapped', help='path to file for unmapped entries',
                        required=False, default="unmapped.gff", type=str)
    parser.add_argument('--useEmergencyMapper', help='use relative mapping with entries that fail to align',
                        required=False, default=True, type=bool)
    parser.add_argument('--specific_entry', help='only translate a specific entry in gff, index (integer). Used for jobs sent to cluster.',
                        required=False, default=None, type=int)
    args = parser.parse_args()
    # important paths
    # ref = {
    #     "fa" : ('/g/steinmetz/project/IESY/genomes/annotations/scramble/genomes/S288C_reference_genome_R64-2-1_20150113/'
    #           'S288C_reference_sequence_R64-2-1_20150113.fsa'),
    #     "gff" : ('/g/steinmetz/project/IESY/genomes/annotations/scramble/genomes/S288C_reference_genome_R64-2-1_20150113/'
    #           'saccharomyces_cerevisiae_R64-2-1_20150113_annotation_only.gff')}
    #
    # query = {
    #     "fa" : ('/g/steinmetz/project/IESY/genomes/annotations/scramble/genomes/'
    #             'JS710.fa'),
    #     "gff" : ('/g/steinmetz/project/IESY/genomes/annotations/scramble/genomes/'
    #             'test_JS710.gff')
    # }
    if args.ref_gff is not None:
        gff = readGFF(args.ref_gff)
        if args.annotation_filter is not None:
            gff = gff[~gff['feature'].isin(args.annotation_filter)]
        if args.specific_entry is None:
            out = mapMultipleGFF(args.ref_fa, gff, args.query_fa, usemp=True, nproc=args.cores,
                use = args.use, blatPath = args.blatpath, unmapped = args.unmapped,
                useEmergencyMapper = args.useEmergencyMapper)
        else:
            if args.specific_entry < gff.shape[0]:
                out = mapGFF(gff_in = gff.iloc[[args.specific_entry]],
                            from_fa = args.ref_fa, to_fa = args.query_fa,
                            use = args.use, blatPath = args.blatpath,
                            unmapped = args.unmapped,
                            useEmergencyMapper = args.useEmergencyMapper, full_gff = gff,
                            onlyHighest=True, correct_tRNA=True, enforceLen=False)
            else:
                print("No such specific entry. Entry outside index.")
                out = None
    else:
        out = None

    if args.segment_table is not None:
        out = addSeg(args.query_fa, args.segment_table, gff = out, usemp=True, nproc=args.cores,
                    use = args.use, blatPath = args.blatpath, unmapped = args.unmapped)
    #out = out.drop_duplicates(["seqname","start","end"])
    if out is not None:
        out = out[~out.index.duplicated(keep='first')]
        #out = out.sort_values(["seqname","start","end"])
        writeGFF(out,args.output_gff)

    # EXAMPLE
    # BASE="/g/steinmetz/project/IESY/genomes/annotations/scramble"
    # REF_FA=$BASE/genomes/S288C_reference_genome_R64-2-1_20150113/S288C_reference_sequence_R64-2-1_20150113.fsa
    # REF_GFF=$BASE/genomes/S288C_reference_genome_R64-2-1_20150113/saccharomyces_cerevisiae_R64-2-1_20150113.gff
    #
    # QUERY_FA_1=/g/steinmetz/project/IESY/genomes/synIXR_scramble/JS734_ERCC92.fasta
    # OUTPUT_GFF_1=$BASE/gff/JS734.gff
    # /g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scripts/scrambleMapper.py \
    # --ref_fa $REF_FA --ref_gff $REF_GFF --query_fa $QUERY_FA_1 --output_gff $OUTPUT_GFF_1
    #
    # QUERY_FA_2=/g/steinmetz/project/IESY/genomes/synIXR_scramble/JS94_ERCC92.fasta
    # OUTPUT_GFF_2=$BASE/gff/JS94.gff
    # /g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scripts/scrambleMapper.py \
    # --ref_fa $REF_FA --ref_gff $REF_GFF --query_fa $QUERY_FA_2 --output_gff $OUTPUT_GFF_2
    #
    # QUERY_FA_3=/g/steinmetz/project/IESY/genomes/synIXR_scramble/JS731_ERCC92.fasta
    # OUTPUT_GFF_3=$BASE/gff/JS710.gff
    # /g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scripts/scrambleMapper.py \
    # --ref_fa $REF_FA --ref_gff $REF_GFF --query_fa $QUERY_FA_3 --output_gff $OUTPUT_GFF_3
    #
    # QUERY_FA_4=/g/steinmetz/project/IESY/genomes/synIXR_scramble/JS731_ERCC92.fasta
    # OUTPUT_GFF_4=$BASE/gff/JS731.gff
    # /g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scripts/scrambleMapper.py \
    # --ref_fa $REF_FA --ref_gff $REF_GFF --query_fa $QUERY_FA_4 --output_gff $OUTPUT_GFF_4
    #
    # QUERY_FA_5=/g/steinmetz/project/IESY/genomes/synIXR_scramble/S288C_ERCC92.fasta
    # OUTPUT_GFF_5=$BASE/gff/S288C.gff
    # /g/steinmetz/brooks/git/steinmetz-lab/yeast2_0/scripts/scrambleMapper.py \
    # --ref_fa $REF_FA --ref_gff $REF_GFF --query_fa $QUERY_FA_5 --output_gff $OUTPUT_GFF_5
