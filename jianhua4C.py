import os
import re
from subprocess import call
from textwrap import dedent
import mappy as mp
import pandas as pd
import sys


def confirm_CWD():
    global cwd
    cwd = os.getcwd()
    if not os.path.exists(f'{cwd}/4Cseqpipe/4cseqpipe.pl'):
        print('wrong working directory')
        exit(1)


def write_4Cseqconf(species):
    if species == 'Mus_musculus':
        track_root = 'trackdb/mm9_trackdb'
    elif species == 'Homo_sapiens':
        track_root = 'trackdb/hg19_trackdb'
    else:
        print('Invalid Species. Use Homo_sapiens or Mus_musculus')
        exit(1)

    conf = f'''\
            index=index.txt
            trackdb_root={track_root}
            trackset_name=4C
            binsize=16
            rawdir=rawdata
            Rscript_path=Rscript'''
    with open(f'{cwd}/4Cseqpipe/4cseqpipe.conf', 'w') as f:
        f.write(dedent(conf))


def write_4Cseqindex(sample_name, species, used_prim, fst_p, sec_p, fst_seq, sec_seq, vp_chr, vp_bp):
    all_index = pd.read_csv(f'{cwd}/4Cseqpipe/index.txt', sep='\t')
    if sample_name in all_index['exp'].values:
        i = all_index[all_index['exp'] == sample_name].iloc[0, 0]
    else:
        i = all_index['id'].max()+1
        if len(used_prim) < 20:
            used_prim = f'{used_prim[0]*(20-len(used_prim))}{used_prim}'
        else:
            used_prim = used_prim[-20:]
        index = {'id': i,
                 'run': sample_name,
                 'lane_no': 0,
                 'exp': sample_name,
                 'primer_seq': used_prim,
                 'species_name': species,
                 'first_cutter_name': fst_p,
                 'first_cutter_seq': fst_seq,
                 'sec_cutter_name': sec_p,
                 'sec_cutter_seq': sec_seq,
                 'linearization_name': None,
                 'linearization_seq': None,
                 'bait_chromo': vp_chr,
                 'bait_coord': vp_bp,
                 'seq_len': 36,
                 'raw_fname': f'{sample_name}.txt'}
        all_index.loc[i] = index
        all_index.fillna('NA').to_csv(
            f'{cwd}/4Cseqpipe/index.txt', sep='\t', index=False)


def run_4Cseq(sample_name, fastq, species, fst_seq, sec_seq):
    if species == 'Mus_musculus':
        track_root = 'trackdb/mm9_trackdb'
    elif species == 'Homo_sapiens':
        track_root = 'trackdb/hg19_trackdb'
    else:
        print('Invalid Species. Use Homo_sapiens or Mus_musculus')
        exit(1)
    write_4Cseqconf(species)
    if not os.path.exists(f'{cwd}/4Cseqpipe/{track_root}/tracks/re/{fst_seq}_dist2{sec_seq}_3prime'):
        os.chdir(f'{cwd}/4Cseqpipe')
        call(f'''perl 4cseqpipe.pl \
            -build_re_db -first_cutter {fst_seq} \
            -second_cutters {sec_seq} \
            -trackdb_root {track_root}''', shell=True)
        os.chdir(cwd)

    os.chdir(f'{cwd}/4Cseqpipe')
    all_index = pd.read_csv(f'{cwd}/4Cseqpipe/index.txt', sep='\t')
    sample_id = all_index[all_index['exp'] == sample_name].iloc[0, 0]
    call(f'''perl 4cseqpipe.pl \
    -fastq2raw \
    -ids {sample_id} \
    -fastq_fn {fastq} \
    -convert_qual 1''', shell=True)

    call(f'''perl 4cseqpipe.pl -map -ids {sample_id}''', shell=True)

    os.chdir(cwd)


def digest_read(name, seq, qual, cutter1, cutter2):
    if cutter1 in seq or cutter2 in seq:
        digest_seq = []
        frg_be = 0
        ith = 0
        for res_enz in re.finditer("|".join([cutter1, cutter2]), seq):
            if ith == 0:
                ith += 1
                frg_en = res_enz.end()
                frg_be = res_enz.start()
                continue
            frg_en = res_enz.end()
            if len(seq[frg_be:frg_en]) >= 20:
                digest_seq.append(
                    f"@{name}\n{seq[frg_be:frg_en]}\n+\n{qual[frg_be:frg_en]}\n")
            frg_be = res_enz.start()
            ith += 1
        if len(seq[frg_be:]) >= 20:
            digest_seq.append(f"@{name}\n{seq[frg_be:]}\n+\n{qual[frg_be:]}\n")
        return '\n'.join(digest_seq)
    else:
        return ''


def dist_fq(sample_name, fst_seq, sec_seq, used_prim, fq1, fq2):
    call(f'mkdir -p {cwd}/dist_fq', shell=True)
    print(f'digest fastq of {sample_name}')
    pipe_fq = open(f'{cwd}/dist_fq/{sample_name}.4cseq.fq', 'w')
    bwa_fq = open(f'{cwd}/dist_fq/{sample_name}.bwa.fq', 'w')
    offset = len(used_prim)-20
    if len(used_prim) < 20:
        used_prim = f'{used_prim[0]*(20-len(used_prim))}{used_prim}'

    for fq in [fq1, fq2]:
        if pd.isnull(fq):
            continue
        for read in mp.fastx_read(fq):
            if read[1].startswith(used_prim):
                if offset >= 0:
                    rd = read[1][offset:offset+36]
                    qual = read[2][offset:offset+36]
                else:
                    rd = ('N'*abs(offset)+read[1])[offset:offset+36]
                    qual = ('F'*abs(offset)+read[1])[offset:offset+36]
                bwa_read = digest_read(
                    read[0], read[1], read[2], fst_seq, sec_seq)
                pipe_read = f'''\
                    @{read[0]}
                    {rd}
                    +
                    {qual}
                    '''
                pipe_fq.write(dedent(pipe_read))
                bwa_fq.write(dedent(bwa_read))
    pipe_fq.close()
    bwa_fq.close()


def align(sample_name, species, vp_chr):
    if species == 'Mus_musculus':
        REF = '~/REF/UCSC/mm9/Sequence/BWAIndex/genome.fa'
    elif species == 'Homo_sapiens':
        REF = '~/REF/UCSC/hg19/Sequence/BWAIndex/genome.fa'
    else:
        print('Invalid Species. Use Homo_sapiens or Mus_musculus')
        exit(1)

    samtools = '~/software/samtools-1.13/samtools'
    call(f'''bwa aln -t 10 \
        {REF} {cwd}/dist_fq/{sample_name}.bwa.fq \
        > {cwd}/alignment/{sample_name}.sai''', shell=True)
    call(f'''bwa samse \
        {REF} {cwd}/alignment/{sample_name}.sai \
        {cwd}/dist_fq/{sample_name}.bwa.fq \
        > {cwd}/alignment/{sample_name}.sam''', shell=True)
    call(f'''{samtools} sort -@ 10 \
    -m 4G -O bam \
    -o {cwd}/alignment/{sample_name}_sorted.bam \
    {cwd}/alignment/{sample_name}.sam''', shell=True)
    call(f'{samtools} index {cwd}/alignment/{sample_name}_sorted.bam', shell=True)

    chrom = f"chr{vp_chr}"
    call(f'''{samtools} view \
    {cwd}/alignment/{sample_name}_sorted.bam {chrom} -O bam \
    -o {cwd}/alignment/{sample_name}_sorted_{chrom}.bam''', shell=True)
    call(f'{samtools} index {cwd}/alignment/{sample_name}_sorted_{chrom}.bam', shell=True)


def plot(sample_name, vp_chr, vp_bp, species, fst_seq, sec_seq, window=100000, ):

    write_4Cseqconf(species)
    os.chdir(f'{cwd}/4Cseqpipe')
    all_index = pd.read_csv(f'{cwd}/4Cseqpipe/index.txt', sep='\t')
    sample_id = all_index[all_index['exp'] == sample_name].iloc[0, 0]
    call(f'''perl 4cseqpipe.pl \
    -nearcis \
    -calc_from {vp_bp-window} \
    -calc_to {vp_bp+window} \
    -stat_type median \
    -trend_resolution 5000 \
    -ids {sample_id} \
    -figure_fn {cwd}/plots/4cpipe/{sample_name}_{window}.pdf''', shell=True)
    os.chdir(cwd)

    if species == 'Mus_musculus':
        track_root = 'trackdb/mm9_trackdb'
    elif species == 'Homo_sapiens':
        track_root = 'trackdb/hg19_trackdb'
    else:
        print('Invalid Species. Use Homo_sapiens or Mus_musculus')
        exit(1)
    if not os.path.exists(f'{cwd}/4Cseqpipe/{track_root}/{fst_seq}_{sec_seq}.txt'):
        call(f'''~/anaconda3/envs/4C/bin/Rscript \
            ./Basic4Cseq_digest.R {fst_seq} {sec_seq} \
            {cwd}/4Cseqpipe/{track_root}/{fst_seq}_{sec_seq}.txt {species}''', shell=True)

    rmap = pd.read_csv(
        f'{cwd}/4Cseqpipe/{track_root}/{fst_seq}_{sec_seq}.txt', sep='\t')
    rmap = rmap[rmap['chromosomeName'] == f'chr{vp_chr}']
    rmap = rmap.iloc[:, :3].copy()
    rmap = rmap.astype({'fragmentStart': int, 'fragmentEnd': int})
    rmap.to_csv(f'{cwd}/coverage/bait.bed',
                sep='\t', index=False, header=False)
    call(
        f'''bedtools coverage \
            -a {cwd}/coverage/bait.bed \
            -b {cwd}/alignment/{sample_name}_sorted.bam \
            > {cwd}/coverage/{sample_name}_cov.txt''', shell=True)
    cov = pd.read_csv(
        f'{cwd}/coverage/{sample_name}_cov.txt', sep='\t', header=None)
    cov[7] = cov[3]/cov[3].sum()*1e6
    cov[[1, 7]].to_csv(
        f'{cwd}/coverage/{sample_name}.txt', sep='\t', index=False, header=False)
    call(f'''~/anaconda3/envs/4C/bin/Rscript \
        {cwd}/peakC.R \
        {cwd}/coverage \
        {cwd}/plots/peakC \
        {sample_name}.txt \
        {vp_bp} {window}''', shell=True)

    call(
        f'echo "chr{vp_chr}	{vp_bp}	{vp_bp}	VP	black" > {cwd}/coverage/vp.bed', shell=True)
    call(f'''~/anaconda3/envs/4C/bin/Rscript Basic4Cseq.R \
        {cwd}/4Cseqpipe/{track_root}/{fst_seq}_{sec_seq}.txt \
        {cwd}/coverage/vp.bed {cwd}/alignment/{sample_name}_sorted_chr{vp_chr}.bam \
        chr{vp_chr} {vp_bp} {window} {sample_name}''', shell=True)


def output_file():
    call(f'echo')


def clean(output_dir):
    for species in ['hg19', 'mm9']:
        call(
            f'rm -rf {cwd}/4Cseqpipe/trackdb/{species}_trackdb/tracks/4C', shell=True)
        call(
            f'mkdir -p {cwd}/4Cseqpipe/trackdb/{species}_trackdb/tracks/4C', shell=True)
    for pipedir in ['figures', 'rawdata', 'stats', 'tables']:
        call(f'rm -rf {cwd}/4Cseqpipe/{pipedir}', shell=True)
        call(f'mkdir -p {cwd}/4Cseqpipe/{pipedir}', shell=True)
    call(f'mkdir -p {output_dir}/alignment', shell=True)
    call(f'mv {cwd}/alignment/*_sorted_chr*.bam* {output_dir}/alignment', shell=True)
    call(f'mkdir -p {output_dir}/coverage', shell=True)
    call(f'mv {cwd}/coverage/*_cov.txt {output_dir}/coverage', shell=True)
    call(f'mv {cwd}/plots {output_dir}/plots', shell=True)
    for workdir in ['alignment', 'coverage', 'dist_fq', 'plots/4cpipe', 'plots/basic4c', 'plots/peakC']:
        call(f'rm -rf {cwd}/{workdir}', shell=True)
        call(f'mkdir -p {cwd}/{workdir}', shell=True)


if __name__ == '__main__':
    confirm_CWD()
    meta = pd.read_excel(sys.argv[1])
    for i in meta.index:
        sample_name, species, fst_p, fst_seq, sec_p, sec_seq, used_prim, unused_prim, vp_chr, vp_bp, fq1, fq2 = meta.iloc[
            i].values
        dist_fq(sample_name, fst_seq, sec_seq, used_prim, fq1, fq2)
        write_4Cseqindex(sample_name, species, used_prim, fst_p,
                         sec_p, fst_seq, sec_seq, vp_chr, vp_bp)
        run_4Cseq(sample_name,
                  f'{cwd}/dist_fq/{sample_name}.4cseq.fq', species, fst_seq, sec_seq)
        align(sample_name, species, vp_chr)
        plot(sample_name, vp_chr, vp_bp, species,
             fst_seq, sec_seq, window=500000,)
    clean(sys.argv[2])
