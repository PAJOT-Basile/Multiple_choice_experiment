import pandas as pd
import os


def get_list_individuals():
    all_barcodes = pd.read_csv(
        "/shared/projects/sexisol/script/Basile/Demultiplexage/barcodes/barcodes_all.tsv", sep="\t").assign(
            libs_index=lambda x: x["ID_Banque"] + "_" + x["code_index"],
    )
    i = 0
    sample_id = ["" for _ in range(len(all_barcodes))]
    for _, row in all_barcodes.iterrows():
        sample_id[i] = "_".join(row.ID_DNA_RAD.split("_")[0:-1])
        i += 1

    all_indivs = all_barcodes["ID_DNA_RAD"]
    all_samples = sample_id
    all_libs_indexes = list(set(all_barcodes["libs_index"]))
    return (all_indivs, all_samples, all_libs_indexes)


def get_projects_lanes():
    indir = "/shared/projects/sexisol/archive/25_ddRAD_novogene/"
    list_projects_lanes = []
    for j in os.listdir(indir):
        if "Jaera" in j:
            for f in os.listdir(indir + j + "/"):
                if "Jaera" in f:
                    list_projects_lanes.append(f)

    return (list_projects_lanes)


def get_lane(wildcards, all_files):
    lanes = [file.split("_")[-1]
             for file in all_files if wildcards.indiv in file]
    return (lanes)
