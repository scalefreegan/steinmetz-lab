#!/usr/bin/env python3

def readFASTA(x, splitKey = None):
    """
    Is sequence file? Load from file if so. File should be FASTA format
    Use pyfasta
    """
    import os
    from pyfaidx import Fasta
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
    import os
    from Bio import SeqIO
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
    import pyfaidx
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
    import pandas as pd
    import numpy as np
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
    import pandas as pd
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
    import os
    import Bio
    import warnings
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

def runMummer(ref, query):
    """
    Align sequences with MUMmer
    """
    from pymummer import coords_file, alignment, nucmer
    import os
    import Bio
    import tempfile
    import datetime
    import pandas as pd
    import pyfaidx
    import uuid

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

def gffMod(df, mummerhits, doReps=False, repNum=1):
    import pandas as pd
    import re
    import urllib
    def add2id(matchobj,i):
        x = str(matchobj.group(0).rstrip(";"))+ "-" + str(i) + ";"
        return x
    out = pd.DataFrame.copy(df)
    seqname = mummerhits["NAME_Q"]
    start_s = int(mummerhits["S2"])
    end_s = int(mummerhits["E2"])
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
        strand = "-"
    else:
        start = start_s
        end = end_s
        strand = "+"
    out.index = pd.MultiIndex.from_tuples([(seqname,start,end)], names=("seqname","start","end"))
    out.loc[(seqname,start,end),"strand"] = strand
    out.loc[(seqname,start,end),"attribute"] = attribute
    return out

def mapGFF(from_fa, gff_in, to_fa, onlyHighest=True, correct_tRNA=True, enforceLen=True):
    import pandas as pd
    from_seq = readFASTA_SeqIO(from_fa)
    feature_seq = getSeq(from_seq,gff_in.index.get_level_values('seqname')[0],
                       gff_in.index.get_level_values('start')[0],
                       gff_in.index.get_level_values('end')[0])
    if feature_seq is not None:
        to_seq = readFASTA_SeqIO(to_fa)
        mummerhits = runMummer(feature_seq,to_seq)
        mummerhits["%_IDY"] = pd.to_numeric(mummerhits["%_IDY"])
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
                if mummerhits.shape[0] > 1:
                    out = gffMod(gff_in, mummerhits.iloc[0], doReps=True, repNum=1)
                    for i in range(1,(mummerhits.shape[0])):
                        #print(i)
                        newrow = gffMod(gff_in,mummerhits.loc[i], doReps=True, repNum=i+1)
                        out = out.append(newrow)
                else:
            	    out = gffMod(gff_in, mummerhits.iloc[0])
                #return mummerhits
                return out
            else:
                return None
        else:
            return None
    else:
        return None

def mapMultipleGFF(from_fa, gff_in, to_fa, usemp=False, nproc=32):
    from tqdm import tqdm
    import multiprocessing as mp
    import pandas as pd
    from itertools import repeat
    if gff_in.shape[0] > 0:
        if usemp:
            print("Using multiprocessing")
            pool = mp.Pool(processes=nproc)
            #results = [pool.apply(mapGFF, args=(from_fa,gff_in.iloc[[x]], to_fa)) for x in range(0,gff_in.shape[0])]
            gff_arg = [gff_in.iloc[[x]] for x in range(0,gff_in.shape[0])]
            results = pool.starmap(mapGFF, zip(repeat(from_fa),gff_arg,repeat(to_fa)))
            #results = pool.map(partial(mapGFF, gff_in=gff_arg), to_fa=to_fa, from_fa=from_fa)
            #results.close()
            #results.join()
            out = pd.concat(results)
            return out
        else:
            print("Using tqdm")
            out = mapGFF(from_fa, gff_in.iloc[[0]], to_fa)
            if gff_in.shape[0] > 1:
                for i in tqdm(range(1,(gff_in.shape[0]))):
                    #print(i)
                    newrow = mapGFF(from_fa, gff_in.iloc[[i]], to_fa)
                    out = out.append(newrow)
            return out
    else:
        return None

if __name__ == "__main__":
    import argparse
    from importlib.machinery import SourceFileLoader
    annotations2keep = [
            'gene',
            'ncRNA_gene',
            'tRNA_gene',
            'centromere',
            'rRNA_gene'
            ]
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_fa', help='file path to reference genome, fasta format',
                        required=True, metavar = "REFERENCE FASTA FILE PATH")
    parser.add_argument('--ref_gff', help='file path to reference genome, gff format',
                        required=True, metavar = "REFERENCE GFF FILE PATH")
    parser.add_argument('--query_fa', help='file path to query genome, fasta format',
                        required=True, metavar = "QUERY FASTA FILE PATH")
    parser.add_argument('--output_gff', help='file path to output gff',
                        required=True, metavar = "OUTPUT GFF FILE PATH")
    parser.add_argument('--annotation_filter', help='list of gff "type" features to keep',
                        required=False, default=annotations2keep, type=list, metavar = "FEATURE TYPE LIST")
    parser.add_argument('--cores', help='number of cores to use',
                        required=False, default=32, type=int)
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
    gff = readGFF(args.ref_gff)
    if len(args.annotation_filter) > 0:
        gff = gff[gff['feature'].isin(args.annotation_filter)]
    out = mapMultipleGFF(args.ref_fa, gff, args.query_fa, usemp=True, nproc=32)
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
