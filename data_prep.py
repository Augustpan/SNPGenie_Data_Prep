#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
import pandas as pd
import os, os.path
import re

def varscan2vcf(fname, output_dir=""):
    '''
        Convert VarScan2 mpileup2snp output format to VCF 4 format

        The converted VCF file contains only AD and DP in INFO and FORMAT columns, 
        which is ready for use of SNPGenie within group analysis.
    '''
    print(fname)
    df = pd.read_csv(fname, sep="\t")

    # VCF headers/metadata
    vcf = [
        '##fileformat=VCFv4.1',
        '##source=VarScan2',
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="">',
        '##INFO=<ID=AD,Number=2,Type=Integer,Description="">',
        '##FILTER=<ID=str10,Description="">',
        '##FILTER=<ID=indelError,Description="">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="">',
        '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(fname)
    ]

    # VCF rows
    for _, row in df.iterrows():
        sp = row[4].split(":")
        fmt = {
            "CHROM": row.Chrom,
            "POS": row.Position,
            "ID": ".",
            "REF": row.Ref,
            "ALT": row.Var,
            "QUAL": ".",
            "FILTER": row[5].split(":")[0],
            "INFO": "DP={};AD={},{}".format(sp[1], sp[2], sp[3]),
            "FORMAT": "DP:AD",
            "SAMPLE": "{}:{},{}".format(sp[1], sp[2], sp[3])
        }
        
        s = "{CHROM}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}".format(**fmt)
        vcf.append(s)

    # Output to VCF file
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    a, _ = os.path.splitext(fname)
    out_fname = os.path.join(output_dir, "{}.vcf".format(a))

    with open(out_fname, "w") as f:
        f.write("\n".join(vcf))

def update_varscan_vcf(fname, output_dir=""):
    with open(fname) as f:
        lines = f.readlines()

    new_lines = []
    for line in lines:
        line = line.strip()
        m = re.match(".+\t\d+\t\.\t\w+\t\w+\t\.\t\w+\t.+\t(.+)\t(.+)", line)
        if m:
            col = m.group(1).split(":")
            val = m.group(2).split(":")
            DP_index = col.index("DP")
            AD_index = col.index("AD")
            
            new_AD_value = "{},{}".format(int(val[DP_index])-int(val[AD_index]), val[AD_index])
            val[AD_index] = new_AD_value
            new_line = re.sub("(.+\t\d+\t\.\t\w+\t\w+\t\.\t\w+\t.+\t.+\t).+", r"\g<1>{}".format(":".join(val)), line)
            new_lines.append(new_line)
        elif line.startswith("#"):
            new_lines.append(line)
        else:
            print("UNRECOGNIZED LINE: ", line)
    
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    with open(os.path.join(output_dir, os.path.basename(fname)), "w") as f:
        f.write("\n".join(new_lines))

def generate_gtf_and_reference(feature_tab="full_feature_tab.csv", input_dir="snpgenie_input", manifest="manifest_vcf", output_dir="snpgenie_input"):
    df = pd.read_csv(feature_tab) 
    df = df[df['4'] != "gene"] # drop all records annotated as gene, only CDS are reserved
    
    if manifest:
        manifest_vcf = {}
        with open(manifest) as f:
            for line in f.readlines():
                a, b = line.strip().split("@")
                manifest_vcf[a] = b
    
    for fname in os.listdir(input_dir): # transverse the input_dir, find all *.vcf files
        if not fname.endswith(".vcf"):
            continue

        with open(os.path.join(input_dir, fname)) as f: # read *.vcf file
            lines = f.readlines()
        
        for line in lines:
            m = re.match("##reference=file://(.+)", line) # find the line indicating location of reference sequence fasta file
            if m:
                ref_path = m.group(1)    
            if not line.startswith("#"):
                seq_name = line.split()[0]
                break
        if manifest:
            m = re.match("VCF(\d+)_SEQ\d+.vcf", os.path.basename(fname))
            ref_path = manifest_vcf[m.group(1)]
        
        slist = []
        for _, row in df[df['1'] == seq_name].iterrows(): # convert feature table to GTF format
            s = "\t".join(
                [fname, 
                 "Python", 
                 "CDS", 
                 str(row['2']), 
                 str(row['3']), 
                 ".", 
                 "+" if row['3'] > row['2'] else "-", 
                 "0", 
                 'gene_id "{}";'.format(row['5'].replace("(","").replace(")",""))])
            slist.append(s)
        
        a, _ = os.path.splitext(fname)
        with open(os.path.join(output_dir, "{}.gtf".format(a)), "w") as f:
            f.write("\n".join(slist))
        
        if not os.path.exists(ref_path): # find reference locally
            ref_path = os.path.join("references", os.path.basename(ref_path))

        for seq in SeqIO.parse(ref_path, "fasta"):
            if seq.id == seq_name:
                SeqIO.write(seq, os.path.join(output_dir, "{}.fasta".format(a)), "fasta")
                break

def split_vcf_varscan(input_dir="varscan_vcfs", output_dir="snpgenie_input"):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".vcf"):
            continue
        
        in_fname = os.path.join(input_dir, fname)
        prefix, _ = os.path.splitext(os.path.basename(fname))

        with open(in_fname) as f:
            lines = f.readlines()
        
        meta = []
        vcf = {}
        for line in lines:
            if line.startswith("#"):
                meta.append(line)
            else:
                seq_name = line.split()[0]
                if seq_name in vcf:
                    vcf[seq_name].append(line)
                else:
                    vcf[seq_name] = [line]
        
        idx = 0
        for key in vcf:
            idx += 1
            text = "".join(meta+vcf[key])
            if vcf[key][-1].startswith("#"):
                print("NO VARIANT WAS CALLED IN: ", key)
                continue
            
            out_fname = os.path.join(output_dir, "VCF{}_SEQ{}.vcf".format(prefix,idx))
            with open(out_fname, "w") as f:
                f.write(text)

def split_vcf_bcftools(input_dir="bcftools_vcfs", output_dir="snpgenie_input"):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for fname in os.listdir(input_dir):
        if not fname.endswith(".vcf"):
            continue
        
        in_fname = os.path.join(input_dir, fname)
        prefix, _ = os.path.splitext(os.path.basename(fname))

        with open(in_fname) as f:
            lines = f.readlines()

        meta = []
        ref = []
        for i, line in enumerate(lines):
            if line.startswith("#"):
                m = re.match("##contig=<ID=([^,]+),length=\d+>", line)
                if m:
                    refseq_id = m.group(1)
                    ref.append((line, refseq_id))
                else:
                    if line.startswith("##reference"):
                        ind_refseq = i
                    meta.append(line)

        vcf = {}
        for line, refseq_id in ref:
            meta_copy = meta.copy()
            meta_copy.insert(ind_refseq+1, line)
            vcf[refseq_id] = meta_copy
        for line in lines:
            if not line.startswith("#"):
                sp = line.split()
                vcf[sp[0]].append(line)

        idx = 0
        for key in vcf:
            idx += 1
            text = "".join(vcf[key])
            if vcf[key][-1].startswith("#"):
                print("NO VARIANT WAS CALLED IN: ", key)
                continue
            
            out_fname = os.path.join(output_dir, "VCF{}_SEQ{}.vcf".format(prefix,idx))
            with open(out_fname, "w") as f:
                f.write(text)

def parse_feature_table(ano):
    df = []
    cache = []
    feature_index = -1
    flag = 1
    for line in ano.split("\n"):
        if line.startswith(">Feature"):
            feature_index += 1 
        else:
            sp = [x.replace(">", "").replace("<", "") for x in line.split()]
            if sp:
                if sp[0] == "product":
                    cache.append("_".join(sp[1:]))
                    flag = 1
                elif sp[0] == "gene":
                    cache.append("_".join(sp[1:]))
                    flag = 1
                elif sp[0] == "note":
                    continue
                else:
                    if flag == 0:
                        cache.append("frame_shift")
                        df.append(cache)
                    cache = [feature_index]
                    cache.extend(sp)
                    flag = 0
                if flag:
                    if cache[-2] not in ["gene", "CDS"]:
                        cache.insert(-1, "CDS_fs")
                    df.append(cache)
    return df

def summarize_feature_table(input_dir = "annotation_files", output="full_feature_tab.csv"):
    FEA = {}
    SEQ = {}
    for fname in os.listdir(input_dir):
        if fname.startswith("Feature"):
            ftname = fname.replace("Feature table file-", "").replace(".txt", "")
            FEA[ftname] = os.path.join(input_dir, fname)

        if fname.startswith("sequence"):
            vname = fname.replace("sequence_", "").replace(".txt", "")
            SEQ[vname] = os.path.join(input_dir, fname)

    df_list = []
    for key in FEA:
        fea_file = FEA[key]
        seq_file = SEQ[key]

        with open(fea_file) as f:
            ano= f.read()

        df = parse_feature_table(ano)

        seq_map = {}
        seq_index = 0
        for seq in SeqIO.parse(seq_file, "fasta"):
            seq_name = seq.description
            seq_map[seq_index] = seq_name
            seq_index += 1

        # in place !!!
        for row in df:
            try:
                seq_map[row[0]]
            except:
                print(row[0])
            row[0] = seq_map[row[0]]
            row.insert(0, key)
            row[1] = row[1].replace(" ", "_")
        df_list.append(pd.DataFrame(df))
    feature_tab = pd.concat(df_list)
    feature_tab.to_csv(output, index=False)
    
    return feature_tab

def rename_fasta_seqs(input_dir="annotation_files", output_dir="references"):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    
    seqs = {}
    for fname in os.listdir(input_dir):
        if fname.startswith("sequence"):
            for seq in SeqIO.parse(os.path.join(input_dir, fname), "fasta"):
                new_name = seq.description.replace(" ", "_")
                seq.description = new_name
                seq.name = new_name
                seq.id = new_name
                
                m = re.match(".+\]_?([a-zA-Z]+_[a-zA-Z]+_[a-zA-Z]+)_.+", seq.id)
                if m:
                    host = m.group(1)
                    if host in seqs:
                        seqs[host].append(seq)
                    else:
                        seqs[host] = [seq]

    for key in seqs:
        SeqIO.write(seqs[key], os.path.join(output_dir, "ref_{}.fa".format(key)), "fasta")

if __name__ == "__main__":
    rename_fasta_seqs()
    summarize_feature_table()

    vc_method = "varscan_vcf"
    if vc_method == "varscan_var":
        var_files = "varscan_vars"
        for fname in os.listdir(var_files):
            if fname.endswith(".var"):
                varscan2vcf(os.path.join(var_files, fname), "varscan_vcfs")
        split_vcf_varscan()
        generate_gtf_and_reference()
    elif vc_method == "varscan_vcf":
        vcf_files = "varscan_vcfs"
        for fname in os.listdir(vcf_files):
            if fname.endswith(".vcf"):
                update_varscan_vcf(os.path.join(vcf_files, fname), output_dir="varscan_vcfs")
        split_vcf_varscan()
        generate_gtf_and_reference()
    elif vc_method == "bcftools":
        split_vcf_bcftools()
        generate_gtf_and_reference(manifest="")
    else:
        print("Invalid vc_method")