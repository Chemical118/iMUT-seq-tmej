import os
import pandas as pd
import numpy as np
import re

samples = ['53bp1_delt', 'artemis_delt', 'blm_delt', 'brca2_delt', 'ctrl_delt', 'exo1_delt', 'fanca_delt', 'lig4_delt', 'pnkp_delt', 'pole_delt', 'poll_delt', 'rad51_delt', 'rad52_delt', 'xrcc4_delt', 'ku70_delt', 'mre11_delt', 'brca1_delt', 'atmrd_delt', 'dnapkd_delt', 'atm_delt', 'atr_delt', 'pold1_delt']
samplesctrl = ['53bp1_ctrl', 'artemis_ctrl', 'blm_ctrl', 'brca2_ctrl', 'exo1_ctrl', 'fanca_ctrl', 'lig4_ctrl', 'pnkp_ctrl', 'pole_ctrl', 'poll_ctrl', 'rad51_ctrl', 'rad52_ctrl', 'xrcc4_ctrl', 'ku70_ctrl', 'mre11_ctrl', 'brca1_ctrl', 'atmrd_ctrl', 'dnapkd_ctrl', 'atm_ctrl', 'atr_ctrl', 'pold1_ctrl']





############################################################################# Read individual replicate mprofiles #############################################################################
############################################################################# Read individual replicate mprofiles #############################################################################
############################################################################# Read individual replicate mprofiles #############################################################################
############################################################################# Read individual replicate mprofiles #############################################################################
############################################################################# Read individual replicate mprofiles #############################################################################

n1 = pd.read_csv("proc/combined_delt1.mprofile", sep="\t")
n2 = pd.read_csv("proc/combined_delt2.mprofile", sep="\t")
n3 = pd.read_csv("proc/combined_delt3.mprofile", sep="\t")
n4 = pd.read_csv("proc/combined_delt4.mprofile", sep="\t")
n5 = pd.read_csv("proc/combined_delt5.mprofile", sep="\t")
n6 = pd.read_csv("proc/combined_delt6.mprofile", sep="\t")
n7 = pd.read_csv("proc/combined_delt7.mprofile", sep="\t")
n8 = pd.read_csv("proc/combined_delt8.mprofile", sep="\t")
n9 = pd.read_csv("proc/combined_delt9.mprofile", sep="\t")
n10 = pd.read_csv("proc/combined_delt10.mprofile", sep="\t")

n1c = pd.read_csv("proc/combined_ctrl1.mprofile", sep="\t")
n2c = pd.read_csv("proc/combined_ctrl2.mprofile", sep="\t")
n3c = pd.read_csv("proc/combined_ctrl3.mprofile", sep="\t")
n4c = pd.read_csv("proc/combined_ctrl4.mprofile", sep="\t")
n5c = pd.read_csv("proc/combined_ctrl5.mprofile", sep="\t")
n6c = pd.read_csv("proc/combined_ctrl6.mprofile", sep="\t")
n7c = pd.read_csv("proc/combined_ctrl7.mprofile", sep="\t")
n8c = pd.read_csv("proc/combined_ctrl8.mprofile", sep="\t")
n9c = pd.read_csv("proc/combined_ctrl9.mprofile", sep="\t")
n10c = pd.read_csv("proc/combined_ctrl10.mprofile", sep="\t")





############################################################################# Mutation averages for output #############################################################################
############################################################################# Mutation averages for output #############################################################################
############################################################################# Mutation averages for output #############################################################################
############################################################################# Mutation averages for output #############################################################################
############################################################################# Mutation averages for output #############################################################################

# This takes the df lists created before and calculates per replicate mutation averages as well as MMEJ rate and deletion lengths
##### MMEJ from delta function #####
# End of the indel sequence should match the sequence after the deletion e.g. -8atgctct would have tctactgcat downstream
def mmejcalc(inputdf):
    outl = list()
    for row in inputdf.itertuples():
        commons = [c for c in row[24].split(",") if c != ""]
        mmrate = 0
        mmejlens = ""
        mmejdellens = ""
        totallen = ""
        if commons:
            indydf = inputdf[row[0]+1:row[0]+300]
            downstream = "".join(indydf[indydf["Chromosome"] == row[1]]["Ref.Base"])
            downl = len(downstream)
            for common in commons:
                indel = common.split(":")[0]
                if indel.startswith("-"):
                    rate = common.split(":")[1]
                    match = re.search(r'\d+', indel)
                    totallen+=(match.group()+",")
                    if row[3] <= 2 and (row[3] + int(match.group())) >= (-2):
                        seq = indel[match.end():]
                        seql = len(seq)
                        if seql <= downl: # if the deletion length is shorter than the remaining downstream nucleotides 
                            mmlen = 0 # start the microhomology length count at 0
                            for j in range(0,seql-1): # for j in the length of the deletion
                                if downl != j+seql: # if the deleted nucleotide position isn't beyond the last downstream nucleotide of the amplicon
                                    if seq[j] == downstream[j+seql]: # if next deleted nucleotide is the same as the next updatream nucleotide 
                                        mmlen+=1 # the potential microhomology increases by 1bp 
                                    else:
                                        break # if the next deleted does not equal the next downstream then the homologous stretch has ended 
                                else:
                                    break
                        else:
                            break

                        if mmlen >= 1:
                            mmejlens+=(str(mmlen)+",")
                            mmejdellens+=(match.group()+",")
                            mmrate+=float(rate)
        

        newrow = list(row)[1:]+[mmejdellens, mmejlens, mmrate, totallen]
        outl.append(newrow)
    mmejdf = pd.DataFrame.from_records(outl, columns=list(inputdf.columns)+["MMEJ.Del.Lengths", "MMEJ.Lengths", "MMEJ.Rate", "Deletion.Lengths"])
    return(mmejdf)

def sigcalc(inputdf):
    sigdic = {"C>A":"C>A", "G>T":"C>A", "C>T":"C>T", "G>A":"C>T", "C>G":"C>G", "G>C":"C>G", "T>A":"T>A", "A>T":"T>A", "T>C":"T>C", "A>G":"T>C", "T>G":"T>G", "A>C":"T>G"}
    outl = list()
    for row in inputdf.itertuples():
        ratdic = {"C>A":0, "C>T":0, "C>G":0, "T>A":0, "T>C":0, "T>G":0}
        mutlis = ["C>A", "C>T", "C>G", "T>A", "T>C", "T>G"]
        bas = row[4]
        arat = row.toA
        trat = row.toT
        grat = row.toG
        crat = row.toC
        for mut, rat in zip(["A", "T", "G", "C"], [arat, trat, grat, crat]):
            if mut != bas:
                ratdic[sigdic[bas+">"+mut]] += rat
        
        newrow = list(row)[1:]+[ratdic[mut] for mut in mutlis]
        outl.append(newrow)

    sigdf = pd.DataFrame.from_records(outl, columns=list(inputdf.columns)+mutlis)
    return(sigdf)

samples = ['53bp1_delt', 'artemis_delt', 'blm_delt', 'brca2_delt', 'ctrl_delt', 'exo1_delt', 'fanca_delt', 'lig4_delt', 'pnkp_delt', 'pole_delt', 'poll_delt', 'rad51_delt', 'rad52_delt', 'xrcc4_delt', 'ku70_delt', 'mre11_delt', 'brca1_delt', 'atmrd_delt', 'dnapkd_delt', 'atm_delt', 'atr_delt', 'pold1_delt']
samplesctrl = ['53bp1_ctrl', 'artemis_ctrl', 'blm_ctrl', 'brca2_ctrl', 'exo1_ctrl', 'fanca_ctrl', 'lig4_ctrl', 'pnkp_ctrl', 'pole_ctrl', 'poll_ctrl', 'rad51_ctrl', 'rad52_ctrl', 'xrcc4_ctrl', 'ku70_ctrl', 'mre11_ctrl', 'brca1_ctrl', 'atmrd_ctrl', 'dnapkd_ctrl', 'atm_ctrl', 'atr_ctrl', 'pold1_ctrl']
samples = ['ctrl_delt', 'ku70_delt', 'dnapkd_delt', 'artemis_delt', '53bp1_delt', 'pnkp_delt', 'poll_delt', 'xrcc4_delt', 'lig4_delt', 'mre11_delt', 'atmrd_delt', 'brca1_delt', 'blm_delt', 'exo1_delt', 'brca2_delt', 'fanca_delt', 'rad52_delt', 'rad51_delt', 'pold1_delt', 'pole_delt', 'atm_delt', 'atr_delt']
samplesctrl = ['ku70_ctrl', 'dnapkd_ctrl', 'artemis_ctrl', '53bp1_ctrl', 'pnkp_ctrl', 'poll_ctrl', 'xrcc4_ctrl', 'lig4_ctrl', 'mre11_ctrl', 'atmrd_ctrl', 'brca1_ctrl', 'blm_ctrl', 'exo1_ctrl', 'brca2_ctrl', 'fanca_ctrl', 'rad52_ctrl', 'rad51_ctrl', 'pold1_ctrl', 'pole_ctrl', 'atm_ctrl', 'atr_ctrl']
samplesraw = ['ctrl', 'dnapkd', 'ku70', '53bp1', 'artemis', 'pnkp', 'poll', 'xrcc4', 'lig4', 'mre11', 'atmrd', 'brca1', 'blm', 'exo1', 'brca2', 'rad51', 'rad52', 'pole', 'fanca', 'atm', 'atr', 'pold1']
allmutations = ['fromA', 'fromT', 'fromG', 'fromC', 'toA', 'toT', 'toG', 'toC', "C>A", "C>T", "C>G", "T>A", "T>C", "T>G", 'Transitions', 'Transversions', 'Total.SNVs', 'Break.Adjacent.SNVs', 'Distant.SNVs', 'Insertions', 'Deletions', 'Small.Indels', 'Mid.Indels', 'Large.Indels', 'Break.Adjacent.Indels', 'Distant.Indels', 'Break.Adjacent.Total', 'Distant.Total', 'Break.Small', 'Break.Large', 'Distant.Small', 'Distant.Large', 'In.Frame', 'Frame.Shift', 'MMEJ.Rate']

##### Build a total reps df #####
dfs = [n2, n3, n4, n5, n6, n7, n8, n9, n10] 
dfcs = [n2c, n3c, n4c, n5c, n6c, n7c, n8c, n9c, n10c] 
for i in range(0,9):
    dfs[i]["Rep"] = i+2
    dfcs[i]["Rep"] = i+2
totaldf = pd.concat(dfs, ignore_index=True)
totaldfc = pd.concat(dfcs, ignore_index=True)





mmejdf = mmejcalc(totaldf)
mmejdfc = mmejcalc(totaldfc)

mmejdf.to_csv("proc/combined_delt_n110_mmej.mprofile", sep="\t", index=False)
mmejdfc.to_csv("proc/combined_ctrl_n110_mmej.mprofile", sep="\t", index=False)

# mmejdf = pd.read_csv("proc/combined_delt_n110_mmej.mprofile", sep="\t")
# mmejdfc = pd.read_csv("proc/combined_ctrl_n110_mmej.mprofile", sep="\t")

sigdf = sigcalc(mmejdf)
sigdfc = sigcalc(mmejdfc)

sigdf.to_csv("proc/combined_delt_n110_mmejsig.mprofile", sep="\t", index=False)
sigdfc.to_csv("proc/combined_ctrl_n110_mmejsig.mprofile", sep="\t", index=False)

# sigdf = pd.read_csv("proc/combined_delt_n110_mmejsig.mprofile", sep="\t")
# sigdfc = pd.read_csv("proc/combined_ctrl_n110_mmejsig.mprofile", sep="\t")










##### Avg mutation rates #####
# Delt
pathmutations = list()
for mut in allmutations:
    pathmutations.append(mut)
    pathmutations.append(mut)
    pathmutations.append(mut)
    pathmutations.append(mut)

outdf = pd.DataFrame(index=pathmutations)
for sample in samples:
    sampdf = sigdf[(sigdf["Sample"]==sample)]
    reps = sorted(list(set(sampdf["Rep"])))
    for rep in range(2,11):
        if rep in reps:
            repdf = sampdf[sampdf["Rep"]==rep]
            mutout = list()
            for mutation in allmutations:
                for path, le in zip(["CTRL", "HR", "NH", ("HR", "NH")], [5, 10, 10, 20]):
                    pathdf = repdf[repdf["Chromosome"].str.startswith(path)]
                    out = sum(pathdf[mutation])/le
                    mutout.append(out)
            outdf[str(sample+str(rep))] = mutout
        else:
            outdf[str(sample+str(rep))] = np.nan
outdf.to_csv("pdout/avg_allmutations_delt_total.csv", index=True)

# Ctrl
outdf = pd.DataFrame(index=pathmutations)
for sample in samplesctrl:
    sampdf = sigdfc[(sigdfc["Sample"]==sample)]
    reps = sorted(list(set(sampdf["Rep"])))
    for rep in range(2,11):
        if rep in reps:
            repdf = sampdf[sampdf["Rep"]==rep]
            mutout = list()
            for mutation in allmutations:
                for path,le in zip(["CTRL", "HR", "NH", ("HR", "NH")], [5, 10, 10, 20]):
                    pathdf = repdf[repdf["Chromosome"].str.startswith(path)]
                    out = sum(pathdf[mutation])/le
                    mutout.append(out)
            outdf[str(sample+str(rep))] = mutout
        else:
            outdf[str(sample+str(rep))] = np.nan
outdf.to_csv("pdout/avg_allmutations_ctrl_total.csv", index=True)





##### Deletion length tables #####
# Total Deletions 
outdf = pd.DataFrame(index=range(0,100000))
for p in ["CT", "HR", "NH", ("HR", "NH")]:
    for sample in samples:
        sampdf = sigdf[(sigdf["Sample"]==sample) & (sigdf["Chromosome"].str.startswith(p))]
        reps = sorted(list(set(sampdf["Rep"])))
        for rep in range(2,11):
            if rep in reps:
                repdf = sampdf[sampdf["Rep"]==rep]
                mutout = list()
                for v in list(repdf.dropna()["Deletion.Lengths"]):
                    if v != "":
                        vs = [int(v3) for v3 in v.split(",") if v3 != "" and int(v3) > 1]
                        for v2 in vs:
                            mutout.append(v2)
                outdf[str(sample+str(rep)+"_"+p[0])] = pd.Series(mutout)
            else:
                outdf[str(sample+str(rep)+"_"+p[0])] = np.nan
outdf.to_csv("pdout/dellens_total.csv", index=False)

outdf = pd.DataFrame(index=range(0,100000))
for p in [("HR", "NH"), "HR", "NH", "CT"]:
    for sample in samples:
        sampdf = sigdf[(sigdf["Sample"]==sample) & (sigdf["Chromosome"].str.startswith(p)) & (sigdf["Rep"] < 8)]
        mutout = list()
        for v in list(sampdf.dropna()["Deletion.Lengths"]):
            if v != "":
                vs = [int(v3) for v3 in v.split(",") if v3 != "" and int(v3) > 1]
                for v2 in vs:
                    mutout.append(v2)
        outdf[str(sample+"_"+p[0]+p[1])] = pd.Series(mutout)

outdf.to_csv("pdout/dellens_total.csv", index=False)


# MMEJ Deletions
outdf = pd.DataFrame(index=range(0,100000))
for p in [("HR", "NH"), "HR", "NH"]:
    for sample in samples:
        sampdf = sigdf[(sigdf["Sample"]==sample) & (sigdf["Chromosome"].str.startswith(p)) & (sigdf["Rep"] < 8)]
        mutout = list()
        for v in list(sampdf.dropna()["MMEJ.Del.Lengths"]):
            if v != "":
                vs = [int(v3) for v3 in v.split(",") if v3 != ""]
                for v2 in vs:
                    mutout.append(v2)
        outdf[str(sample+"_"+p[0]+p[1])] = pd.Series(mutout)

outdf.to_csv("pdout/mmej_dellens_total.csv", index=False)

# MMEJ homologies
outdf = pd.DataFrame(index=range(0,100000))
for p in [("HR", "NH"), "HR", "NH"]:
    for sample in samples:
        sampdf = sigdf[(sigdf["Sample"]==sample) & (sigdf["Chromosome"].str.startswith(p)) & (sigdf["Rep"] < 8)]
        mutout = list()
        for v in list(sampdf.dropna()["MMEJ.Lengths"]):
            if v != "":
                vs = [int(v3) for v3 in v.split(",") if v3 != ""]
                for v2 in vs:
                    mutout.append(v2)
        outdf[str(sample+"_"+p[0]+p[1])] = pd.Series(mutout)

outdf.to_csv("pdout/mmej_homlens_total.csv", index=False)





##### Site Averages #####
# Delt
sitemuts = ["Total.SNVs", 'Break.Adjacent.SNVs', 'Distant.SNVs', 'Insertions', 'Deletions', 'Small.Indels', 'Mid.Indels', 'Large.Indels', 'Break.Adjacent.Indels', 'Distant.Indels', 'Break.Adjacent.Total', 'Distant.Total', 'Break.Small', 'Break.Large', 'Distant.Small', 'Distant.Large', "MMEJ.Rate"]
outmuts = ["totsnv", 'adjsnv', 'distsnv', 'ins', 'del', 'sml', 'mid', 'lrg', 'adjind', 'distind', 'adjtot', 'disttot', 'adjsml', 'adjlrg', 'distsml', 'distlrg', "mmej"]
for mutation, mutstr in zip(sitemuts, outmuts):
    sites=sorted(list(set(sigdf["Chromosome"])))
    outdf = pd.DataFrame(index=sites)
    for sample in samples:
        sampdf = sigdf[sigdf["Sample"]==sample]
        reps = sorted(list(set(sampdf["Rep"])))
        for rep in range(2,11):
            if rep in reps:
                repdf = sampdf[sampdf["Rep"]==rep]
                mutout = list()
                for site in sites:
                    sitedf = repdf[repdf["Chromosome"] == site]
                    out = sum(sitedf[mutation])
                    mutout.append(out)
                outdf[str(sample+str(rep))] = mutout
            else:
                outdf[str(sample+str(rep))] = np.nan
    outdf.to_csv("pdout/site_avgs/site_average_delt_"+mutstr+".csv", index=True)

    # Ctrl
    sites=sorted(list(set(sigdfc["Chromosome"])))
    outdf = pd.DataFrame(index=sites)
    for sample in samplesctrl:
        sampdf = sigdfc[sigdfc["Sample"]==sample]
        reps = sorted(list(set(sampdf["Rep"])))
        for rep in range(2,11):
            if rep in reps:
                repdf = sampdf[sampdf["Rep"]==rep]
                mutout = list()
                for site in sites:
                    sitedf = repdf[repdf["Chromosome"] == site]
                    out = sum(sitedf[mutation])
                    mutout.append(out)
                outdf[str(sample+str(rep))] = mutout
            else:
                outdf[str(sample+str(rep))] = np.nan
    outdf.to_csv("pdout/site_avgs/site_average_ctrl_"+mutstr+".csv", index=True)


