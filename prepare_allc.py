#!/usr/bin/env python3
import os
import gzip
import math
import time
import click
import shlex
import subprocess

import numpy as np
import pandas as pd

from scipy.stats import mannwhitneyu

WORKDIR = "/baimoc/wangjing/ngs/meth/pipeline/results/11_methylpy"
GENOMEFILE = "/baimoc/wangjing/genome/bovine/bosTau9.fa"

FILES = {
    "adult_forebrain": {
        "Abrain_1": 0.0062,
        "Abrain_2": 0.0059,
    },
    "adult_muscle": {
        "Ajirou_1": 0.0024,
        "Ajirou_2": 0.0035,
    },
    "adult_rumen": {
        "Awei_1": 0.0017,
        "Awei_2": 0.0045,
    },
    "fetal_rumen": {
        "wei": 0.0034,
    },
    "fetal_heart": {
        "xin_1": 0.00285,
    },
    "fetal_kidney": {
        "shen_1": 0.00286,
    },
    "fetal_liver": {
        "gan_1": 0.00292,
    },
    "fetal_lung": {
        "fei_1": 0.00304,
    },
    "fetal_testis": {
        "gaowan": 0.0031,
    },
    "fetal_forebrain": {"qiannao_1": 0.0031, "qiannao_2": 0.0032},
    "fetal_hindbrain": {"hounao_1": 0.0031, "hounao_2": 0.0032},
    "fetal_muscle": {
        "jirou_1": 0.0032,
        "jirou_2": 0.0031,
    },
    "bESCs_F7": {"F7_P11": 0.0033, "F7_5_22": 0.0031, "F7_50_16": 0.0029},
    # "bESCs_F7_5": {
    #     "F7_5_22": 0.0031,
    # },
    # "bESCs_F7_50": {
    #     "F7_50_16": 0.0029,
    # },
    "bEPSCs_B18": {
        "B18_p32": 0.00276,
        "B18_P35": 0.00284,
    },
    "bEPSCs_AGS": {
        "AGS_bEPSC_P40": 0.00289,
        "AGS_bEPSC_p39": 0.0032,
    },
    "adult_Adipose": {
        "Adipose_1": 0.0021,
        "Adipose_2": 0.0036,
    },
    "adult_Blood": {
        "Blood_1": 0.0035,
        "Blood_2": 0.0036,
    },
    "adult_Heart": {
        "Heart_1": 0.0026,
        # "Heart_2": 0,
    },
    "adult_Kidney": {
        "Kidney_1": 0.0044,
        # "Kidney_2": 0,
    },
    "adult_Lung": {
        "Lung_1": 0.0049,
        "Lung_2": 0.0019,
    },
    "adult_Mammary": {
        "Mammary_1": 0.0057,
        "Mammary_2": 0.0059,
    },
    "adult_Spleen": {
        "Spleen_1": 0.0038,
    },
    "adult_Liver": {
        "Liver_1": 0.007,
    },
    "adult_Uterus": {"Uterus_1": 0.0069},
    "adult_Ileum": {"Ileum_1": 0.0033},
}


def open_allc_file(allc_file: str):
    if allc_file.endswith("gz"):  # so it works for .gz and .bgz
        f = gzip.open(allc_file, "rt")
    else:
        f = open(allc_file, "r")
    return f


def test_allc(sample: str, rate: float):
    path = os.path.join(WORKDIR, "allc")
    file = os.path.join(path, sample + ".methylpy.tsv.gz")
    if rate > 0:
        cmd = f"methylpy test-allc --allc-file {file} --sample {sample} --unmethylated-control {rate} --num-procs 8 --min-cov 3 --remove-chr-prefix True --path-to-output {path}"
        print(cmd)
        print("==")
        subprocess.check_call(shlex.split(cmd))
    else:
        outfile = os.path.join(path, "allc_" + sample + ".tsv.gz")
        cmd = f"cp {file} {outfile}"
        subprocess.check_call(shlex.split(cmd))


def filter_methylKit(sample: str):
    f = open_allc_file(
        os.path.join(WORKDIR, "cytosine_report", sample + ".cytosine_report.txt.gz")
    )
    g = gzip.open(
        os.path.join(WORKDIR, "cytosine_report", sample + ".methylkit.gz"), "wt"
    )
    for line in f:
        fields = line.rstrip().split("\t")
        # print(fields)
        if fields[5] == "CG":
            g.write("\t".join(fields) + "\n")


def format_allc(sample: str):
    f = open_allc_file(
        os.path.join(WORKDIR, "cytosine_report", sample + ".cytosine_report.txt.gz")
    )
    g = gzip.open(os.path.join(WORKDIR, "allc", sample + ".methylpy.tsv.gz"), "wt")
    for line in f:
        fields = line.rstrip().split("\t")
        g.write(
            fields[0][3:]
            + "\t"
            + fields[1]
            + "\t"
            + fields[2]
            + "\t"
            + fields[6]
            + "\t"
            + fields[3]
            + "\t"
            + str(int(fields[3]) + int(fields[4]))
            + "\t1"
            + "\n"
        )
    f.close()
    g.close()


def filter_CG(sample: str):
    infile = os.path.join(WORKDIR, "allc", "allc_" + sample + ".tsv.gz")
    outfile = os.path.join(WORKDIR, "allc", "allCG", "allCG_" + sample + ".tsv.gz")
    cmd = f"methylpy filter-allc --allc-file {infile} --output-file {outfile} --mc-type CGN --min-cov 3 --compress-output True --num-procs 8"
    print(cmd)
    print("==")
    subprocess.check_call(shlex.split(cmd))


def filter_CH(sample: str):
    infile = os.path.join(WORKDIR, "allc", "allc_" + sample + ".tsv.gz")
    outfile = os.path.join(WORKDIR, "allc", "allCH", "allCH_" + sample + ".tsv.gz")
    cmd = f"methylpy filter-allc --allc-file {infile} --output-file {outfile} --mc-type CHN --min-cov 3 --compress-output True --num-procs 8"
    print(cmd)
    print("==")
    subprocess.check_call(shlex.split(cmd))


def allc_to_bigwig(sample: str):

    infile = os.path.join(WORKDIR, "allc", "allc_" + sample + ".tsv.gz")
    outfile1 = os.path.join(WORKDIR, "bigwig", "allCG" + sample + ".bw")
    outfile2 = os.path.join(WORKDIR, "bigwig", "allCH" + sample + ".bw")

    # for s in FILES.keys():
    #     if sample in FILES[s].keys():
    #         nonCov_rate = FILES[s][sample]

    cmd1 = f"methylpy allc-to-bigwig --allc-file {infile} --output-file {outfile1} --ref-fasta {GENOMEFILE} --mc-type CGN --add-chr-prefix True"
    print(cmd1)
    subprocess.check_call(shlex.split(cmd1))

    cmd2 = f"methylpy allc-to-bigwig --allc-file {infile} --output-file {outfile2} --ref-fasta {GENOMEFILE} --mc-type CHN --add-chr-prefix True"
    print(cmd2)
    subprocess.check_call(shlex.split(cmd2))


def merge_allc(name: str, samples: list):
    if len(samples) > 1:
        infiles = " ".join(
            [
                os.path.join(WORKDIR, "allc", "allCG", f"allCG_{s}.tsv.gz")
                for s in samples
            ]
        )
        outfile = os.path.join(
            WORKDIR, "allc", "allCG", "merged", f"allCG_{name}.tsv.gz"
        )
        cmd = f"methylpy merge-allc --allc-files {infiles} --output-file {outfile} --num-procs 8 --compress-output True"
    else:
        infile = os.path.join(WORKDIR, "allc", "allCG", f"allCG_{samples[0]}.tsv.gz")
        outfile = os.path.join(
            WORKDIR, "allc", "allCG", "merged", f"allCG_{name}.tsv.gz"
        )
        cmd = f"cp {infile} {outfile}"

    if not os.path.exists(outfile):
        print(cmd)
        print("==")
        subprocess.check_call(shlex.split(cmd))


def merge_allch(name: str, samples: list):
    if len(samples) > 1:
        infiles = " ".join(
            [
                os.path.join(WORKDIR, "allc", "allCH", f"allCH_{s}.tsv.gz")
                for s in samples
            ]
        )
        outfile = os.path.join(
            WORKDIR, "allc", "allCH", "merged", f"allCH_{name}.tsv.gz"
        )
        cmd = f"methylpy merge-allc --allc-files {infiles} --output-file {outfile} --num-procs 8 --compress-output True"
    else:
        infile = os.path.join(WORKDIR, "allc", "allCH", f"allCH_{samples[0]}.tsv.gz")
        outfile = os.path.join(
            WORKDIR, "allc", "allCH", "merged", f"allCH_{name}.tsv.gz"
        )
        cmd = f"cp {infile} {outfile}"

    if not os.path.exists(outfile):
        print(cmd)
        print("==")
        subprocess.check_call(shlex.split(cmd))


def allc_DMR(exclude: str, outdir: str = ""):
    files = ""
    samples = ""
    excludes = exclude.split(",")
    for s in FILES.keys():
        merge_allc(s, list(FILES[s].keys()))
        # merge_allch(s, list(FILES[s].keys()))

        if s not in excludes:
            file = os.path.join(WORKDIR, "allc", "allCG", "merged", f"allCG_{s}.tsv.gz")
            files = files + " " + file
            samples = samples + " " + s
    outpath = os.path.join(WORKDIR, "CG_DMR", outdir)
    cmd = f"methylpy DMRfind --allc-files {files.strip()} --samples {samples.strip()} --mc-type CGN --num-procs 8 --output-prefix {outpath} --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 X --min-num-dms 3"
    print(cmd)
    print("==")
    subprocess.check_call(shlex.split(cmd))


def steps_allc(step, sample=None):
    if step == "-1":
        return
    elif step == "0":
        for s in FILES[sample].keys():
            filter_methylKit(s)
            format_allc(s)
            test_allc(s, FILES[sample][s])
            filter_CG(s)
            filter_CH(s)
            allc_to_bigwig(s)
    elif step == "1":
        for s in FILES[sample].keys():
            filter_methylKit(s)
            format_allc(s)
    elif step == "2":
        for s in FILES[sample].keys():
            test_allc(s, FILES[sample][s])
    elif step == "3":
        for s in FILES[sample].keys():
            filter_CG(s)
            filter_CH(s)
    elif step == "4":
        for s in FILES[sample].keys():
            allc_to_bigwig(s)
    else:
        print("ERROR: Invalid step number.")
        exit(1)


@click.group()
def main():
    pass


@main.command()
@click.option(
    "--steps",
    type=click.Choice(["-1", "0", "1", "2", "3", "4"]),
    default="-1",
    help="step to process: -1: do nothing; 0: all; 1: format allc; 2: test allc; 3: filter allc; 4: allc to bigwig.",
)
@click.option(
    "--sample", default="ALL", help="sample to process, [ALL] for all samples."
)
def process(steps, sample):
    if sample in FILES.keys():
        print("Sample: ", sample)
        print(time.asctime(time.localtime(time.time())) + "\n")
        steps_allc(steps, sample)
        print("Done.")
    else:
        if sample == "ALL":
            # if dmr:
            #     print("Perform DMRfind for all samples.")
            #     print(time.asctime(time.localtime(time.time())) + "\n")
            #     allc_DMR(exclude)
            # else:
            for name in FILES.keys():
                print("Sample: ", name)
                print(time.asctime(time.localtime(time.time())) + "\n")
                steps_allc(steps, name)
            print("Done.")
        else:
            print("Invalid sample name.")
            exit(1)


@main.command()
@click.option(
    "--exclude", default="", help="samples to exclude for DMRfind, comma separated."
)
@click.option("--prefix", default="", help="output prefix for DMR result files.")
def cg_dmr(exclude, prefix):
    print(time.asctime(time.localtime(time.time())) + "\n")
    allc_DMR(exclude, prefix)

    outpath = prefix + "_rms_results_collapsed.tsv"
    dmr = pd.read_table(outpath, sep="\t")

    dmr["hypermethylated_samples"] = ""
    dmr["hypomethylated_samples"] = ""
    dmr.insert(6, "baseline", np.nan)
    dmr.head()

    new_dmr = dmr.apply(tissue_specific_DMR, axis=1)
    new_dmr.to_csv(prefix + ".cg_dmr.csv")


def tissue_specific_DMR(cur_line: pd.Series, difference: float = 0.3) -> pd.Series:
    tissue_up = []
    tissue_down = []
    line = cur_line.copy().iloc[7:].sort_values()
    n = math.ceil(len(line) / 2)
    if pd.isna(line).sum() < n:
        idx = (
            line.dropna()
            .rolling(n, min_periods=n)
            .agg(lambda x: x.iloc[-1] - x.iloc[0])
            .dropna()
            .reset_index(drop=True)
            .idxmin()
        )

        baseline = line.iloc[idx : idx + n].mean()

        for key, value in line.items():
            if value - baseline >= difference:
                tissue_up.append(key.replace("methylation_level_", ""))
            elif baseline - value >= difference:
                tissue_down.append(key.replace("methylation_level_", ""))
        cur_line["hypermethylated_samples"] = ",".join(tissue_up)
        cur_line["hypomethylated_samples"] = ",".join(tissue_down)
        cur_line["baseline"] = baseline
    return cur_line


@main.command()
@click.option("--prefix", default="", help="prefix to add.")
@click.option("--suffix", default="", help="suffix to add.")
@click.option("--merge", default=False, help="merge rep or not.")
def print_all_samples(prefix, suffix, merge):
    if merge:
        for s in FILES.keys():
            print(prefix + s + suffix, end=" ")
    else:
        for s in FILES.keys():
            for sample in FILES[s].keys():
                print(prefix + sample + suffix, end=" ")
    print("\n")


def get_samples_mCH(sample1, sample2, mCH_level, outdir):
    outfile1 = mCH_level.iloc[:, :3].copy()
    outfile2 = mCH_level.iloc[:, :3].copy()

    if sample1 in FILES.keys() and sample2 in FILES.keys():
        samples1 = list(FILES[sample1].keys())
        samples2 = list(FILES[sample2].keys())

        if len(samples1) > 1:
            mCH_level_sample1 = []
            for sample in samples1:
                if "methylation_level_" + sample in mCH_level.columns:
                    tmp = (
                        mCH_level["methylation_level_" + sample]
                        - FILES[sample1][sample]
                    )
                    # tmp = mCH_level["methylation_level_" + sample].copy()
                    tmp[tmp < 0] = 0

                    med = tmp.dropna().median()
                    if med > 0:
                        tmp = tmp / med
                    else:
                        raise Exception("median == 0")
                else:
                    raise Exception("Invalid column name.")

                mCH_level_sample1.append(tmp)

            outfile1[sample1] = np.mean(mCH_level_sample1, axis=0)
            for i in range(len(mCH_level_sample1)):
                outfile2[mCH_level_sample1[i].name] = mCH_level_sample1[i]
        else:
            if "methylation_level_" + samples1[0] in mCH_level.columns:
                tmp = (
                    mCH_level["methylation_level_" + samples1[0]]
                    - FILES[sample1][samples1[0]]
                )
                # tmp = mCH_level["methylation_level_" + samples1[0]].copy()
                tmp[tmp < 0] = 0
                med = tmp.dropna().median()
                tmp = tmp / med
                outfile1[sample1] = tmp
                outfile2["methylation_level_" + samples1[0]] = tmp
            else:
                raise Exception("Invalid column name.")
        print("check sample1 done.")

        if len(samples2) > 1:
            mCH_level_sample2 = []
            for sample in samples2:
                if "methylation_level_" + sample in mCH_level.columns:
                    tmp = (
                        mCH_level["methylation_level_" + sample]
                        - FILES[sample2][sample]
                    )
                    # tmp = mCH_level["methylation_level_" + sample].copy()
                    tmp[tmp < 0] = 0

                    med = tmp.dropna().median()
                    if med > 0:
                        tmp = tmp / med
                    else:
                        raise Exception("median == 0")
                else:
                    raise Exception("Invalid column name.")

                mCH_level_sample2.append(tmp)
            outfile1[sample2] = np.mean(mCH_level_sample2, axis=0)
            for i in range(len(mCH_level_sample2)):
                outfile2[mCH_level_sample2[i].name] = mCH_level_sample2[i]
        else:
            if "methylation_level_" + samples2[0] in mCH_level.columns:
                tmp = (
                    mCH_level["methylation_level_" + samples2[0]]
                    - FILES[sample2][samples2[0]]
                )
                # tmp = mCH_level["methylation_level_" + samples2[0]].copy()
                tmp[tmp < 0] = 0
                med = tmp.dropna().median()
                tmp = tmp / med
                outfile1[sample2] = tmp
                outfile2["methylation_level_" + samples2[0]] = tmp
            else:
                raise Exception("Invalid column name.")
        print("check sample2 done.")
    outfile2.to_csv(
        outdir + "/" + sample1 + "_vs_" + sample2 + ".normalized_mCH.bed",
        sep="\t",
        index=False,
    )
    ch2bw(outdir + "/" + sample1 + "_vs_" + sample2 + ".normalized_mCH.bed", outdir)
    return outfile1


def ch2bw(file, outdir):
    df = pd.read_table(file, sep="\t", low_memory=False)

    for x in df.columns[3:]:
        print("Normalized sample to bigwig: ", x.replace("methylation_level_", ""))
        name = outdir + "/" + x.replace("methylation_level_", "") + ".tmp.bedGraph"
        bw_name = outdir + "/" + x.replace("methylation_level_", "") + ".bw"
        tmp_df = df.loc[:, ["chrom", "start", "end", x]].copy().dropna()
        tmp_df["chrom"] = tmp_df["chrom"].apply(lambda x: "{}{}".format("chr", x))
        tmp_df.to_csv(name, sep="\t", index=False, header=False)
        cmd = f"bedGraphToBigWig {name} {GENOMEFILE}.fai {bw_name}"
        subprocess.check_call(shlex.split(cmd))
        subprocess.check_call(["rm", name])


@main.command()
@click.option("--file", help="mCH methylation level file path.")
@click.option("--sample1", help="sample control for mCH mega DMR test.")
@click.option("--sample2", help="sample treat for mCH mega DMR test.")
@click.option("--pvalue", default=0.05, help="p value cut off(fdr).")
@click.option("--fc", default=1.0, help="log fold change cut off(absolute value).")
@click.option("--outdir", default=".", help="output dir.")
def ch_dmr(file, sample1, sample2, pvalue, fc, outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    print(time.asctime(time.localtime(time.time())) + "\n")

    mCH_level = pd.read_table(file, sep="\t", low_memory=False)
    outfile = mCH_level.iloc[:, :3].copy()

    if sample1 in FILES.keys() and sample2 in FILES.keys():
        outfile = get_samples_mCH(sample1, sample2, mCH_level, outdir)

        outfile["p_value"] = np.nan
        outfile["fold_change"] = 0
        n = outfile.shape[0]

        pval_list = []

        print("start mann-whitney u test for 10 adjacent bins.")
        for i in range(0, n, 10):
            sample_inputs_1 = outfile.iloc[i : i + 10, 3]
            sample_inputs_2 = outfile.iloc[i : i + 10, 4]

            if (
                sample_inputs_1.isnull().sum() <= 4
                and sample_inputs_2.isnull().sum() <= 4
            ):
                res = mannwhitneyu(
                    sample_inputs_1.dropna(),
                    sample_inputs_2.dropna(),
                    use_continuity=True,
                )

                pval_list.append(res.pvalue)
                if n - i >= 10:
                    for j in range(i, i + 10):
                        outfile.iloc[j, 5] = res.pvalue
                        if outfile.iloc[j, 4] >= 0 and outfile.iloc[j, 3] >= 0:
                            outfile.iloc[j, 6] = np.log2(
                                (outfile.iloc[j, 4] + 0.0000001)
                                / (outfile.iloc[j, 3] + 0.0000001)
                            )
                        else:
                            outfile.iloc[j, 6] = 0
                else:
                    for j in range(i, n):
                        outfile.iloc[j, 5] = res.pvalue
                        if outfile.iloc[j, 4] >= 0 and outfile.iloc[j, 3] >= 0:
                            outfile.iloc[j, 6] = np.log2(
                                (outfile.iloc[j, 4] + 0.0000001)
                                / (outfile.iloc[j, 3] + 0.0000001)
                            )
                        else:
                            outfile.iloc[j, 6] = 0
            else:
                pval_list.append(1)
                # outfile.iloc[i, 5] = 1
        print("compute p adjust using fdr.")

        # p_adjust = list(fdr(pval_list))
        # pval_list["p_adjust"] = p_adjust

        # new_file = pd.merge(
        #     outfile, pval_list, on=["chrom", "start", "end", "p_value"], how="left"
        # )
        outfile["p_adjust"] = [val for val in fdr(pval_list) for i in range(10)][:n]
        print("all done.")
        outfile.to_csv(
            outdir + "/" + sample1 + "_vs_" + sample2 + ".ch_dmr.bed",
            sep="\t",
            index=False,
        )

        # up
        up = pd.DataFrame(data=None, columns=outfile.columns)
        # down
        down = pd.DataFrame(data=None, columns=outfile.columns)
        for _, row in outfile.iterrows():
            if float(row["p_adjust"]) <= pvalue:
                if float(row["fold_change"]) > fc:
                    up.loc[len(up)] = row
                elif float(row["fold_change"]) < -fc:
                    down.loc[len(up)] = row

        if up.shape[0] > 0 and down.shape[0] > 0:
            up["chrom"] = ["chr" + str(x) for x in up["chrom"]]
            down["chrom"] = ["chr" + str(x) for x in down["chrom"]]
            up.to_csv(
                outdir + "/" + sample1 + "_vs_" + sample2 + ".up.bed",
                sep="\t",
                index=False,
                header=False,
            )
            down.to_csv(
                outdir + "/" + sample1 + "_vs_" + sample2 + ".down.bed",
                sep="\t",
                index=False,
                header=False,
            )
        else:
            print("no difference, pass.")
    else:
        raise Exception("Invalid sample name.")


def fdr(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


if __name__ == "__main__":
    main()
