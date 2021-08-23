#!/usr/bin/env python
# coding: utf-8

from Bio import SeqIO
import pandas as pd
import os, os.path
import re

def generate_gtf_and_reference(feature_tab="full_feature_tab.csv", input_dir="snpgenie_input", output_dir="snpgenie_input"):
    df = pd.read_csv(feature_tab) 
    df = df[df['4'] != "gene"] # drop all records annotated as gene, only CDS are reserved
    
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
                 'gene_id "{}";'.format(row['5'])])
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

def split_vcf(input_dir="vcf_files", output_dir="snpgenie_input"):
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
    split_vcf()
    generate_gtf_and_reference()