# Without a genetic algorithm these parameters at this granularity would take 1,239,300,000,000 to complete
# With less granularity, good results could be obtained in 6,048,000,000 alignments

from __future__ import division
import subprocess as sb
import sys
from math import log
import numpy as np
import os

try:
    os.mkdir(sys.argv[3])
except OSError:
    pass

def dict_sort(dic):
    return([k for k, v in sorted(dic.items(), key=lambda item: item[1], reverse=True)])

def aligner(dictarg, f1=sys.argv[1], f2=sys.argv[2], sam=sys.argv[3]):
    argue = {key:val for key, val in dictarg.items()}
    defaults = {"l": "32", "dpad": "15", "gbar": "4", "mp": "6,2", "rdg":"5,3", "rfg":"5,3", "scoremin":"-0.6,-0.6"}
    keys = set(argue.keys())
    for key, val in defaults.items():
        if key not in keys:
            argue[key] = val
    proc = sb.Popen("bowtie2 -t -p 60 --fr --maxins 400 --no-discordant --no-mixed -D 100 -R 50 -N 1 --ignore-quals --no-1mm-upfront --np 0 -L "+argue["l"]+" --dpad "+argue["dpad"]+" --gbar "+argue["gbar"]+" --mp "+argue["mp"]+" --rdg "+argue["rdg"]+" --rfg "+argue["rfg"]+" --score-min L,"+argue["scoremin"]+" -x /mnt/data/R11/raw_sequencing_data/Aldo_raw_seq/genomes/mutseq_25_0221/mutseq_25_0221 -1 "+f1+" -2 "+f2+" -S ./"+sam+"/"+sam+"_align.sam", shell=True, stderr=sb.PIPE)
    err = proc.communicate()[1]
    if "exited with value" not in err:
        total = int(err.split("\n")[4].split(" ")[0])
        exact = float(int(err.split("\n")[7].split(" ")[4])/total*100)
        multi = float(int(err.split("\n")[8].split(" ")[4])/total*100)
        non = float(int(err.split("\n")[6].split(" ")[4])/total*100)
        raw_time = err.split("\n")[3].split(" ")[-1]
        mins = (int(raw_time.split(":")[0])*60) + int(raw_time.split(":")[1]) + (int(raw_time.split(":")[2])/60)
        fitness = log(1/(non+0.00001)) - (2*multi)
        return(fitness, [str(exact), str(multi), str(exact+multi), str(non), str(mins)], err)
    else:
        return(-1000, ["0", "0", "0", "0", "0"], err)

def selection(indict, iters=2):
    fitdict = {k:v for k,v in indict.items()}
    tops = list()
    count = 0
    while count < iters:
        val_list = [val[0] for val in fitdict.values()]
        if len(val_list) == 0:
            tops.append(tops[-1])
            count+=1
        else:
            max_val = max(val_list)
            filtered = {k: v for k,v in fitdict.items() if v[0] == max_val}
            if len(filtered) == 1:
                top = list(filtered.keys())[0]
                tops.append(top)
                del(fitdict[top])
                count+=1
            else:
                exacts = dict()
                times = dict()
                for key, val in filtered.items():
                    exacts[key] = float(val[1][0])
                    times[key] = float(val[1][-1])
                max_exact = max(exacts.values())
                min_time = min(times.values())
                filtered2 = {k: v for k,v in filtered.items() if float(v[1][0]) == max_exact}
                if len(filtered2) == 1:
                    top = list(filtered2.keys())[0]
                    tops.append(top)
                    del(fitdict[top])
                    count+=1
                else:
                    filtered3 = {k: v for k,v in filtered2.items() if float(v[1][-1]) == min_time}
                    top = list(filtered3.keys())[0]
                    tops.append(top)
                    del(fitdict[top])
                    count+=1
    return(tuple(tops))

def combolist(l1, l2, greater=False):
    outlist = list()
    for m1 in l1:
        for m2 in l2:
            if greater:
                if m2 <= m1:
                    if (str(m1)+","+str(m2)) not in outlist:
                        outlist.append(str(m1)+","+str(m2))
            else:
                if (str(m1)+","+str(m2)) not in outlist:
                    outlist.append(str(m1)+","+str(m2))

    return(outlist)

def valar(parames, windowdict, top_dict): 
    ranges_func = {}
    for param in parames:
        top_val = top_dict[param]
        if param in ["mp", "rdg", "rfg"]:
            top_split = [float(val) for val in top_val.split(",")]
            window1 = windowdict[str(param +"1")]
            range_temp1 = [str(val) for val in np.arange(top_split[0]-window1[0], top_split[0]+window1[0]+window1[1], window1[1]) if float(val) >= 0]
            window2 = windowdict[str(param +"2")]
            range_temp2 = [str(val) for val in np.arange(top_split[1]-window2[0], top_split[1]+window2[0]+window2[1], window2[1]) if float(val) > 0]
            if param == "mp":
                ranges_func[param] = combolist(range_temp1, range_temp2, greater=True)
            else:
                ranges_func[param] = combolist(range_temp1, range_temp2)
        elif param == "scoremin":
            top_split = [float(val) for val in top_val.split(",")]
            window1 = windowdict[str(param +"1")]
            range_temp1 = [str(val) for val in np.arange(top_split[0]-window1[0], top_split[0]+window1[0]+window1[1], window1[1]) if float(val) <= 0]
            window2 = windowdict[str(param +"2")]
            range_temp2 = [str(val) for val in np.arange(top_split[1]-window2[0], top_split[1]+window2[0]+window2[1], window2[1]) if float(val) <= 0]
            ranges_func[param] = combolist(range_temp1, range_temp2)
        else:
            top_val = float(top_val)
            window = windowdict[param]
            range_temp = [str(val) for val in np.arange(top_val-window[0], top_val+window[0]+window[1], window[1]) if val >= 0]
            ranges_func[param] = range_temp
    return(ranges_func)

total_log = open(sys.argv[3]+"/total_log.txt", "w")
total_log.write("Parmeter\tFitness\tExact\tMultimaps\tTotal\tUnaligned\tTime\n")
total_log.flush()
results_log = open(sys.argv[3]+"/results_log.txt", "w")

combined_log = open("machine1_results.csv", 'a')

print('\n')
i=0
sys.stdout.write("\rAlignments run: "+str(i)+" ")
sys.stdout.flush()


# First get the defualt baseline alignment rate
base_fitness, reports, output = aligner(dict())
total_log.write("Baseline" +"\t"+ str(base_fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
i+=1
results_log.write("Baseline alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
sys.stdout.write("\rAlignments run: "+str(i)+" ")
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()

combined_log.write(sys.argv[3]+","+ str(base_fitness) +","+ ",".join(reports) +",")









params = ["l", "gbar", "mp", "rdg", "rfg"]

ranges = {}
ranges["l"] = ["15", "20", "25", "30", "32"]
ranges["gbar"] = ["1", "2", "4"]
mp1 = ["1", "3"]
mp2 = ["0.5", "1", "1.5", "2"]
ranges["mp"] = combolist(mp1, mp2, greater=True)
rdg1 = ["0", "2", "4"]
rdg2 = ["1", "3"]
ranges["rdg"] = combolist(rdg1, rdg2)
rfg1 = ["0", "2", "4"]
rfg2 = ["1", "3"]
ranges["rfg"] = combolist(rfg1, rfg2)



# First round we take measurements for each paramter at low granularity.
# We need to take the stderr results and extract the alignment frequency as well as the multi-map frequency.
# Build dicts for the top two results for each parameter and pass them to the next round.

r1a_results = {}
r1b_results = {}
for param in params:
    values = ranges[param]
    param_results = {}
    temp_dict = {}
    for val in values:
        temp_dict[param] = val
        fitness, reports, output = aligner(temp_dict)
        if fitness != -1000:
            param_results[val] = [fitness, reports]
        total_log.write(param +"_"+ str(val) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
        total_log.flush()
        i+=1
        sys.stdout.write("\rAlignments run: "+str(i)+" ")
        sys.stdout.flush()
        align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
        align_i_file.write(output)
        align_i_file.close()
    selected = selection(param_results)
    r1a_results[param] = selected[0]
    r1b_results[param] = selected[1]



result_text = [k+"-"+str(v) for k, v in r1a_results.items()]
results_log.write("Round 1 complete in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound one complete.\n")
sys.stdout.flush()










# Now we have the top two values for each parameter.
# We then need to compare these two to each other with the other paramteres set to their best option.
# This is in case parameter changes alter the optimal values of other parameters.

# First get the result with every paramter at top value.
r1_temp = {k:v for k,v in r1a_results.items()}
base_fitness, reports, output = aligner(r1_temp)
total_log.write("Round1" +"\t"+ str(base_fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
i+=1
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
sys.stdout.write("\rAlignments run: "+str(i)+" ")
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()



# Then test changing each parameter to their second value and check if this is better.
# If it is, the top value and second value are switched.
# This will continue to run until no more changes are made.
count = 1
paramchecks = list()
while count > 0:
    count = 0
    r1_temp = {k:v for k,v in r1a_results.items()}
    base_fitness, reports, output = aligner(r1_temp)
    for param in params:
        temp_values = r1_temp
        temp_values[param] = r1b_results[param]
        if temp_values not in paramchecks:
            paramchecks.append(temp_values)
            fitness, reports, output = aligner(temp_values)
            total_log.write(param +"_"+ str(temp_values[param]) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
            total_log.flush()
            i+=1
            sys.stdout.write("\rAlignments run: "+str(i))
            sys.stdout.flush()
            align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
            align_i_file.write(output)
            align_i_file.close()
            if fitness > base_fitness:
                old_top = r1a_results[param]
                r1a_results[param] = r1b_results[param]
                r1b_results[param] = old_top
                count+=1



result_text = [k+"-"+str(v) for k, v in r1a_results.items()]
results_log.write("Round 1 optimised in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound one optimised.\n")
sys.stdout.flush()

# Run the optimised results together to determine optimisation effectiveness.
fitness, reports, output = aligner(r1a_results)
total_log.write("Round1_optimised" +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
i+=1
sys.stdout.write("\rAlignments run: "+str(i))
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()










# Now we likely have the best parameters, however they lack granularity.
# We need to repeat this whole process, but with a range of values for each parameter around the top value.

# First, define the new granular ranges of values.
windows = {}
windows["l"] = [3, 1]
windows["gbar"] = [1, 1]
windows["mp1"] = [0.4, 0.2]
windows["mp2"] = [0.4, 0.2]
windows["rdg1"] = [1, 1]
windows["rdg2"] = [1, 1]
windows["rfg1"] = [1, 1]
windows["rfg2"] = [1, 1]
ranges2 = valar(params, windows, r1a_results)

params = params + ["dpad", "scoremin"]

scoremin1 = ["-20", "-10", "-5", "0"]
scoremin2 = ["-2", "-1", "-0.5", "0"]
ranges2["scoremin"] = combolist(scoremin1, scoremin2)   
ranges2["dpad"] = ["5", "10", "20", "30", "50"]

# Then we can repeat the selection process using the previous optimal results as a baseline.
r2a_results = {}
r2b_results = {}
for param in params:
    old_vals = list()
    count = 0
    while count <= 3:
        values = ranges2[param]
        if values == old_vals:
            count+=4
        else:
            param_results = {}
            temp_dict = {k:v for k,v in r1a_results.items()}
            for val in values: 
                temp_dict[param] = val
                fitness, reports, output = aligner(temp_dict)
                if fitness != -1000:
                    param_results[val] = [fitness, reports]
                total_log.write(param +"_"+ str(val) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
                total_log.flush()
                i+=1
                sys.stdout.write("\rAlignments run: "+str(i)+" ")
                sys.stdout.flush()
                align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
                align_i_file.write(output)
                align_i_file.close()

            selected = selection(param_results)
            if selected[0] == ranges2[param][0] or selected[0] == ranges2[param][-1]:
                count+=1
                r1a_results[param] = selected[0]
                ranges2 = valar(params[:-2], windows, r1a_results)
                ranges2["scoremin"] = combolist(scoremin1, scoremin2)   
                ranges2["dpad"] = ["5", "10", "20", "30", "50"]
                old_vals = values
                total_log.write("re window!\n")
                total_log.flush()
            else:
                count+=4
    r2a_results[param] = selected[0]
    r2b_results[param] = selected[1]



result_text = [k+"-"+str(v) for k, v in r2a_results.items()]
results_log.write("Round 2 complete in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound two complete.\n")
sys.stdout.flush()










# Once again we have the top two values for each parameter.
# We repeat the same process as before to determine if these are optimal when combined  with each other.

# First get the result with every paramter at top value.
r2_temp = {k:v for k,v in r2a_results.items()}
base_fitness, reports, output = aligner(r2_temp)
total_log.write("Round2" +"\t"+ str(base_fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
i+=1
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
sys.stdout.write("\rAlignments run: "+str(i)+" ")
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()



# Then test changing each parameter to their second value and check if this is better.
# If it is, the top value and second value are switched.
# This will continue to run until no more changes are made.

count = 1
paramchecks = list()
while count > 0:
    count = 0
    r2_temp = {k:v for k,v in r2a_results.items()}
    base_fitness, reports, output = aligner(r2_temp)
    for param in params:
        temp_values = {k:v for k,v in r2_temp.items()}
        temp_values[param] = r2b_results[param]
        if temp_values not in paramchecks:
            paramchecks.append(temp_values)
            fitness, reports, output = aligner(temp_values)
            total_log.write(param +"_"+ str(temp_values[param]) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
            total_log.flush()
            i+=1
            sys.stdout.write("\rAlignments run: "+str(i))
            sys.stdout.flush()
            align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
            align_i_file.write(output)
            align_i_file.close()
            if fitness > base_fitness:
                old_top = r2a_results[param]
                r2a_results[param] = r2b_results[param]
                r2b_results[param] = old_top
                count+=1



result_text = [k+"-"+str(v) for k, v in r2a_results.items()]
results_log.write("Round 2 optimised in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound two optimised.\n")
sys.stdout.flush()

# Run the optimised results together to determine optimisation effectiveness.
fitness, reports, output = aligner(r2a_results)
total_log.write("Round2_optimised" +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
i+=1
sys.stdout.write("\rAlignments run: "+str(i))
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()






 



# We're almost there, one last round to increase the granularity 
# We need to repeat this whole process, but with a range of values for each parameter around the top value.

# First, define the new granular ranges of values.
windows = {}
windows["l"] = [3, 1]
windows["gbar"] = [1, 1]
windows["mp1"] = [0.15, 0.05]
windows["mp2"] = [0.1, 0.02]
windows["rdg1"] = [1, 1]
windows["rdg2"] = [1, 1]
windows["rfg1"] = [1, 1]
windows["rfg2"] = [1, 1]
windows["scoremin1"] = [5, 1]
windows["scoremin2"] = [0.5, 0.1]
windows["dpad"] = [5, 1]
ranges3 = valar(params, windows, r2a_results)



# Then we can repeat the selection process using the previous optimal results as a baseline.
r3a_results = {}
r3b_results = {}
r3c_results = {}
for param in params:
    old_vals = list()
    count = 0
    while count <= 3:
        values = ranges3[param]
        if values == old_vals:
            count+=4
        else:
            param_results = {}
            temp_dict = {k:v for k,v in r2a_results.items()}
            for val in values: 
                temp_dict[param] = val
                fitness, reports, output = aligner(temp_dict)
                if fitness != -1000:
                    param_results[val] = [fitness, reports]
                total_log.write(param +"_"+ str(val) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
                total_log.flush()
                i+=1
                sys.stdout.write("\rAlignments run: "+str(i)+" ")
                sys.stdout.flush()
                align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
                align_i_file.write(output)
                align_i_file.close()
            
            selected = selection(param_results, iters=3)
            if selected[0] == ranges3[param][0] or selected[0] == ranges3[param][-1]:
                count+=1
                r2a_results[param] = selected[0]
                ranges3 = valar(params, windows, r2a_results)
                old_vals = values
            else:
                count+=4
    r3a_results[param] = selected[0]
    r3b_results[param] = selected[1]
    r3c_results[param] = selected[2]



result_text = [k+"-"+str(v) for k, v in r3a_results.items()]
results_log.write("Round 3 complete in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound three complete.\n")
sys.stdout.flush()










# Once again we have the top two values for each parameter.
# We repeat the same process as before to determine if these are optimal when combined  with each other.

# First get the result with every paramter at top value.
r3_temp = {k:v for k,v in r3a_results.items()}
base_fitness, reports, output = aligner(r3_temp)
total_log.write("Round3" +"\t"+ str(base_fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
i+=1
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
sys.stdout.write("\rAlignments run: "+str(i)+" ")
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()

# Then test changing each parameter to their second value and check if this is better.
# If it is, the top value and second value are switched.
# This will continue to run until no more changes are made.

count=1
paramchecks = list()
while count >0:
    count=0
    r3_temp = {k:v for k,v in r3a_results.items()}
    base_fitness, reports, output = aligner(r3_temp)
    for param in params:
        temp_values = {k:v for k,v in r3_temp.items()}
        temp_valuesc = {k:v for k,v in r3_temp.items()}
        temp_values[param] = r3b_results[param]
        temp_valuesc[param] = r3c_results[param]
        if temp_values not in paramchecks or temp_valuesc not in paramchecks:
            paramchecks.append(temp_values)
            paramchecks.append(temp_valuesc)
            fitness, reports, output = aligner(temp_values)
            total_log.write(param +"_"+ str(temp_values[param]) +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
            total_log.flush()
            i+=1
            sys.stdout.write("\rAlignments run: "+str(i))
            sys.stdout.flush()
            align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
            align_i_file.write(output)
            align_i_file.close()

            fitnessc, reportsc, outputc = aligner(temp_values)
            total_log.write(param +"_"+ str(temp_values[param]) +"\t"+ str(fitnessc) +"\t"+ "\t".join(reportsc) + "\n")
            total_log.flush()
            i+=1
            sys.stdout.write("\rAlignments run: "+str(i))
            sys.stdout.flush()
            align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
            align_i_file.write(output)
            align_i_file.close()

            if fitnessc > fitness:
                old_top = r3b_results[param]
                r3b_results[param] = r3c_results[param]
                r3c_results[param] = old_top
                fitness = fitnessc
                reports = reportsc
                output = outputc

            if fitness > base_fitness:
                old_top = r3a_results[param]
                r3a_results[param] = r3b_results[param]
                r3b_results[param] = old_top
                count+=1



result_text = [k+"-"+str(v) for k, v in r3a_results.items()]
results_log.write("Round 3 optimised in "+str(i)+" alignments.\nThe optimised parameters are:\n")
results_log.write(", ".join(result_text)+"\n")
results_log.flush()

print("\rRound three optimised.\n")
sys.stdout.flush()

# Run the optimised results together to determine optimisation effectiveness.
fitness, reports, output = aligner(r3a_results)
total_log.write("Round3_optimised" +"\t"+ str(fitness) +"\t"+ "\t".join(reports) + "\n")
total_log.flush()
results_log.write("Combined alignment frequency: " + reports[2] + "\nMultimapping frequency: " + reports[1] + "\nCompletion time: " + reports[-1] + "\n\n\n")
results_log.flush()
i+=1
sys.stdout.write("\rAlignments run: "+str(i))
sys.stdout.flush()
align_i_file = open(sys.argv[3]+"/alignment_"+str(i)+".txt", "w")
align_i_file.write(output)
align_i_file.close()


combined_log.write(sys.argv[3]+","+ str(base_fitness) +","+ ",".join(reports) + "\n")
combined_log.close()

final_param_file = open("machine1_params.txt", 'a')
final_param_file.write(sys.argv[3]+"\t")
for param in params:
    final_param_file.write(str(r3a_results[param])+"\t")
final_param_file.write("\n")
final_param_file.close()






# Time to wrap the files up and end the script...

results_log.close()
total_log.close()
print("\rOptimisation complete.\n")
sys.stdout.flush()
print("\rTotal alignments run: "+str(i)+"\n")
sys.stdout.flush()
