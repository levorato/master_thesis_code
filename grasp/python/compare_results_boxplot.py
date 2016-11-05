# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

# <editor-fold desc="Command-line arguments to calculate chinese box plot graph.">
# /home/mlevorato/Downloads/mestrado-output-temp/output/ils-seq-SetParA/agents-seq-l1-ils10s-a0.4-fi-imb-2h
# /home/mlevorato/Downloads/mestrado-output-temp/output/grasp-bullmpi/agents-seq-l1-grasp400s-a1.0-fi-imb-2h
# --labels
# SeqILS
# SeqGRASP
# --excludefiles
# c4n16k16pin0.7p-0.0p+0.0.g
# c4n32k20pin0.7p-0.6p+0.6.g
# c4n64k32pin0.8p-0.0p+0.1.g
# c4n64k32pin0.8p-0.0p+0.4.g
# c4n64k32pin0.8p-0.0p+0.6.g
# c4n64k32pin0.8p-0.0p+0.8.g
# c4n64k32pin0.8p-0.0p+1.0.g
# c4n64k32pin0.8p-0.2p+1.0.g
# c4n64k32pin0.8p-0.6p+0.0.g
# c4n64k32pin0.8p-0.8p+0.0.g
# c4n64k32pin0.8p-1.0p+0.0.g
# </editor-fold>

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import argparse
import pandas as pd
import scipy.stats as stats
import numpy as np
import random
import matplotlib as mpl
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from matplotlib.backends.backend_pdf import PdfPages

## agg backend is used to create plot as a .png file
mpl.use('agg')

import matplotlib.pyplot as plt

def main(argv):

    csv.field_size_limit(100000)

    parser = argparse.ArgumentParser(description='Compare GRASP / ILS - CC result files, generating box plot graphs.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the result files (one for each experiment / MH)')
    parser.add_argument('--filefilter', default='*.csv', required=False,
                        help='the file extension for result files (default: *.csv)')
    parser.add_argument('--labels', nargs='+',
                        help='the experiment labels (one for each experiment / MH)')
    parser.add_argument('--excludefiles', nargs='+',
                        help='instance file name(s) to exlude from boxplot graphs')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter
    labels = args.labels
    exclude = args.excludefiles
    #sol_multiplier = 1.001
    #sol_multiplier = 1.0016#18123
    sol_multiplier = 1.0

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter
    print 'Instance files to be excluded are ', exclude

    processCCResult(folders, labels, sol_multiplier, exclude)


def processCCResult(folders, labels, sol_multiplier, exclude):

    # for each folder (algorithm), determines the best solution value
    # compares the solution values of each folder and uses the worst solution value as Target I(P)
    # Using this Target I(P), determines the average time spent by each algorithm and outputs a table with the comparison
    # Also outputs a sequence of all execution times of each algorithm, to allow the creation of TTTplots

    bestSolValues = []
    worstSolValues = []

    for folder in folders:
        #print "Processing folder " + ''.join(folder)
        # CC results
        worst_file_summary = dict()
        best_file_summary = dict()
        previous_filename = ""

        for root1, subFolders1, files1 in os.walk(folder):
            # sort dirs and files
            subFolders1.sort()
            files1.sort()

            #print "Processing folder " + ''.join(root1)
            if (''.join(root1) != folder):  # process only subfolders
                current_folder = (root1[root1.rfind(os.path.sep) + 1 :])
                if '.g' in current_folder:
                    filename = current_folder

                    for root, subFolders, files in os.walk(root1):
                        # Process CC results
                        if os.path.isfile(root + os.path.sep + "cc-result.txt"):
                            filename = (root[:root.rfind(os.path.sep)])
                            filename = filename[filename.rfind(os.path.sep)+1:]
                            datetime = root[root.rfind(os.path.sep) + 1:]

                            best_value = float(100000000)

                            content_file = open(root + os.path.sep + "cc-result.txt", 'r')
                            try:
                                content = content_file.read()

                                reader = csv.reader(StringIO.StringIO(content), delimiter='=')
                                for row in reader:
                                    linestring = ''.join(row)
                                    column = []
                                    for col in row:
                                        column.append(col)
                                    if linestring.startswith('I(P)'):
                                        # obtains the best result found by a specific execution of a specific node (can be parallel)
                                        best_value = float(str(column[1]).strip())
                                        break
                            finally:
                                content_file.close()
                            # captura o PIOR resultado dadas todas as execucoes de um mesmo grafo / instancia
                            # deve ser o pior resultado encontrado, dadas todas as execucoes de um determinado algoritmo
                            if best_file_summary.has_key(filename):
                                current_value = float(best_file_summary[filename])
                                if (best_value < current_value):
                                    best_file_summary[filename] = str(best_value)
                            else:
                                best_file_summary[filename] = str(best_value)

                            if worst_file_summary.has_key(filename):
                                current_value = float(worst_file_summary[filename])
                                if (best_value > current_value):
                                    worst_file_summary[filename] = str(best_value)
                            else:
                                worst_file_summary[filename] = str(best_value)
                    # end process all result files of an instance / filename
                # end loop
                # process last file
        print "\nProcessing CC Results...\n"
        bestSolValues.append(best_file_summary)
        worstSolValues.append(worst_file_summary)

        # Save CC results of all executions of all instances to csv file: Time to reach best solution and total time spent
        #result_file = open(folder + os.path.sep + "summary.csv", "w")
        #print "Instance, I(P), I(P)+, I(P)-, k, Iter, Local time(s), Global time(s), Params, Total Iter, Total Comb"
        #result_file.write("ExecutionID; I(P); I(P)+; I(P)-; k; IterBestSol; TimeToBestSol; Global time(s); Params; Total Iter; NumVisitedSolutions\n")
        #for key in sorted(all_files_summary.iterkeys()):
            #print "%s, %s" % (key, all_files_summary[key])
        #    result_file.write("%s; %s\n" % (key, all_files_summary[key]))
        #result_file.close()

    # compare the best results of each algorithm, choosing the worst
    worse_sol = dict()
    best_sol = dict()

    for worst_result_dict in worstSolValues:
        for key in sorted(worst_result_dict.iterkeys()):
            instance_name = str(key)
            sol_value = float(worst_result_dict[key])
            if worse_sol.has_key(instance_name):
                current_value = float(worse_sol[instance_name])
                if sol_value > current_value:
                    worse_sol[instance_name] = str(sol_value)
            else:
                worse_sol[instance_name] = str(sol_value)

    for best_result_dict in bestSolValues:
        for key in sorted(best_result_dict.iterkeys()):
            instance_name = str(key)
            sol_value = float(best_result_dict[key])
            if best_sol.has_key(instance_name):
                current_value = float(best_sol[instance_name])
                if sol_value < current_value:
                    best_sol[instance_name] = str(sol_value)
            else:
                best_sol[instance_name] = str(sol_value)

    print "\n\nBest solution values for each instance, comparing all algorithms:"
    InstanceTargetSolDataSet = list(zip(best_file_summary.iterkeys(), best_file_summary.itervalues()))
    df = pd.DataFrame(data = InstanceTargetSolDataSet, columns=['Instance', 'Target I(P)'])
    print df

    print "\n\nWorst solution values for each instance, comparing all algorithms:"
    InstanceTargetSolDataSet2 = list(zip(worst_file_summary.iterkeys(), worst_file_summary.itervalues()))
    df2 = pd.DataFrame(data = InstanceTargetSolDataSet2, columns=['Instance', 'Target I(P)'])
    print df2

    # for each instance in **worse_sol**, discover the time each algorithm took to reach that solution
    chsol = dict()
    for key in best_sol.iterkeys():
        chsol[key] = (float(worse_sol[key])*(sol_multiplier))
    ttt_for_algorithm = dict()  # contains the ttt dict for a specific algorithm / folder
    for folder in folders:
        ttt_for_algorithm[folder] = processTimeToTarget(folder, chsol)  # worse_sol, best_sol

    # now, for each instance, we need to calculate the columns: AvgTimeToTarget, Stddev(TTT), Speedup.
    # and also the Standard Error of the Mean (SEM)
    AvgTimeToTarget = dict()
    StddevTTT = dict()
    Anova = []
    Tukey = []
    num_exec = dict()

    Speedup = []

    for folder in folders:
        AvgTimeToTarget[folder] = []
        StddevTTT[folder] = []
        num_exec[folder] = []

    for instance_name in best_file_summary.iterkeys():
        # for each algorithm,
        count = 0
        data_to_plot = []
        tttseries = []
        for folder in folders:
            tttdict = dict(ttt_for_algorithm[folder])
            tttlist = tttdict[instance_name]
            tttseries.append(tttlist)
            # obtains a list of execution times (time-to-target) of the algorithm for this specific instance
            # calculates the average, stddev and confidence internal for the ttt variable
            avg_ttt =  np.mean(tttlist)
            AvgTimeToTarget[folder].append(avg_ttt)
            stddev_ttt = (np.std(tttlist, ddof=1)/np.sqrt(len(tttlist)) * 1.96) # numpy.std(tttlist)
            StddevTTT[folder].append(stddev_ttt)
            num_exec[folder].append(len(tttlist))

            # student's t-test for time-to-target
            #(tstat, pvalue) = stats.ttest_1samp(list(tttlist), avg_ttt)
            SEM = stats.sem(list(tttlist))
            data_to_plot.append(list(tttlist))

            # print 'For %s solved by %s : StdErrorOfMean = %6.4f' % (instance_name, labels[count], SEM)
            count += 1

        # To run Wilcoxon signed rank test, sample sizes must be equal!
        if len(tttseries[0]) > len(tttseries[1]):
            while len(tttseries[0]) > len(tttseries[1]):
                tttseries[0].pop(random.randint(0,len(tttseries[0])-1))
        elif len(tttseries[1]) > len(tttseries[0]):
            while len(tttseries[1]) > len(tttseries[0]):
                tttseries[1].pop(random.randint(0,len(tttseries[1])-1))


        # Wilcoxon.append(stats.mannwhitneyu(tttseries[0], tttseries[1]))
        # Applies ANOVA to verify that the samples differ
        f, p = stats.f_oneway(tttseries[0], tttseries[1])
        if p > 0.05:
            print "WARN: p-value above 0.05! Samples may not differ!"
        Anova.append((f, p))

        # runs tukey's test to compare all the samples
        ttt_tukey = []
        for item in tttseries[0]:
            ttt_tukey.append((labels[0], item))
        for item in tttseries[1]:
            ttt_tukey.append((labels[1], item))

        names = ['id','ttt']
        formats = ['U30','f8']
        dtype = dict(names = names, formats=formats)
        #array = np.rec.array(ttt_tukey, dtype=dtype)

        #res = pairwise_tukeyhsd(array['ttt'], array['id'], alpha=0.05)
        #res.plot_simultaneous()
        #plt.show()
        #print res
        #print res.meandiffs[0]
        #Tukey.append(res.meandiffs[0])
        Tukey.append(0)

        # box plot of time-to-target
        # Create a figure instance
        #fig = plt.figure(1, figsize=(9, 6))

        # Create an axes instance
        #ax = fig.add_subplot(111)
        #ax.set_xticklabels(list(labels))

        # Create the boxplot
        #bp = ax.boxplot(list(data_to_plot))
        #plt.show()

        # Save the figure
        #fig.savefig(str(instance_name) + '-boxplot.png', bbox_inches='tight')

    # TODO export TTT series to be used in TTT-plot-latex
    print '\nTTT data for use with TTT-Latex\n'
    data_to_plot = {}
    instance_names = []
    for instance_name in best_file_summary.iterkeys():
        if exclude != None:
            if instance_name in exclude:
                print "Skipping file " + str(instance_name) + ':\n'
                continue
        print str(instance_name) + ':\n'
        folder_count = 0
        data_to_plot[instance_name] = []
        instance_names.append(instance_name)
        for folder in folders:
            tttdict = dict(ttt_for_algorithm[folder])
            tttlist = tttdict[instance_name]
            print str(folder[folder.rfind('/') + 1:]) + ': ' + str(tttlist).replace(',', '').replace('[', '').replace(']', '')
            # NEW: uses the time-to-target series in tttseries[0] and tttseries[1] to generate the box plots
            ## combine these different collections into a list
            print "Aggregating box plots for time-to-target series..."
            # Remove the outliers from the list of TTT
            tttarray = np.array(tttlist)
            filtered = tttarray[~is_outlier(tttarray)]
            data_to_plot[instance_name].append(filtered)
            folder_count += 1

    experiment_name = ''
    for folder in folders:
        experiment_name = experiment_name + folder[folder.rfind('/') + 1:] + '-vs-'
    tdir = './compare_results'
    if not os.path.exists(tdir):
        os.makedirs(tdir)
    result_file_prefix = tdir + '/' + str(experiment_name)

    print "Generating box plot graph for time-to-target series for each instance..."
    generate_box_plot(data_to_plot, instance_names, labels, result_file_prefix)
    print "Box plots successfully generated."

    cols = ['Instance', 'Target I(P)']  # 'AvgTimeToTarget', 'Stddev-TTT', 'Speedup'
    CompleteInstanceDataSet = {'Instance' : list(worst_file_summary.iterkeys()), 'Target I(P)' : list(chsol.itervalues()),
                               'ANOVA' : list(Anova), 'Tukey' : list(Tukey) }
    count = 0
    for folder in folders:
        colname = 'AvgTimeToTarget-' + str(labels[count])
        CompleteInstanceDataSet[colname] = AvgTimeToTarget[folder]
        cols.append(colname)
        colname = 'Stddev-TTT-' + str(labels[count])
        CompleteInstanceDataSet[colname] = StddevTTT[folder]
        cols.append(colname)
        colname = 'NumExec-' + str(labels[count])
        CompleteInstanceDataSet[colname] = num_exec[folder]
        cols.append(colname)
        count += 1

    df = pd.DataFrame(CompleteInstanceDataSet)
    print df[cols]


    df.to_csv(result_file_prefix + '.csv')


def processTimeToTarget(folder, worse_sol):

    print "\nProcessing time-to-target for result folder = " + str(folder) + "\n"
    instance_dict = dict()
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        #print "Processing folder " + ''.join(root)
        if (len(files) and ''.join(root) != folder):
            file_list = []
            file_list.extend(glob.glob(root + os.path.sep + "CC*-iterations.csv"))
            file_list.extend(glob.glob(root + os.path.sep + "Node*-iterations.csv"))
            count = len(file_list) - 1

            # Process CC results
            if os.path.isfile(root + os.path.sep + "cc-result.txt"):
                filename = (root[:root.rfind(os.path.sep)])
                datetime = root[root.rfind(os.path.sep) + 1:]
                filename = filename[filename.rfind(os.path.sep) + 1:]
                myfolder = root[:root.rfind(os.path.sep)]
                #print "Processing instance " + myfolder + ", " + filename
                if instance_dict.has_key(filename):
                    if instance_dict[filename] <> myfolder:
                        print "WARN: problem with result folder structure!\n"
                else:
                    instance_dict[filename] = root[:root.rfind(os.path.sep)]
    # end result folder scan
    tttdict = dict()  # tttdict contains the time-to-target list for each instance for this algorithm
    for key in sorted(instance_dict.iterkeys()):
        instance_name = str(instance_dict[key])
        print "Processing instance " + str(instance_name) + ", " + str(key)
        tttlist = processTimeToTargetOnInstance(instance_name, str(key), worse_sol)
        tttdict[key] = tttlist

    return tttdict

def processTimeToTargetOnInstance(folder, filename, worse_sol):
    # CC results
    all_files_summary = dict()
    full_iteration_data = dict()
    target_iteration_history = dict()

    print "\n\nTime-to-target for instance " + str(filename)
    target_value = float(worse_sol[filename])
    if filename.find("-size8000-part0.g") > 0:
        target_value = float(16068.0)
    print "\nDesired target I(P) = " + str(target_value)

    previous_filename = ""
    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()

        #print "Processing folder " + ''.join(root)
        if (len(files) and ''.join(root) != folder):
            file_list = []
            file_list.extend(glob.glob(root + os.path.sep + "CC*-iterations.csv"))
            file_list.extend(glob.glob(root + os.path.sep + "Node*-iterations.csv"))
            count = len(file_list) - 1

            # Process CC results
            if os.path.isfile(root + os.path.sep + "cc-result.txt"):
                filename = (root[:root.rfind(os.path.sep)])
                datetime = root[root.rfind(os.path.sep) + 1:]
                filename = filename[filename.rfind(os.path.sep) + 1:]
                # print "Processing instance " + filename

                value = 0
                pos_value = 0
                neg_value = 0
                K = 0
                iteration = 0
                time = 0

                while count >= 0:
                    # Process all files from all metaheuristic executions, including those in parallel
                    # print "Processing file " + file_list[count] + "\n"
                    content_file = open(file_list[count], 'r')
                    full_iteration_data.clear()
                    try:
                        content = content_file.read()

                        reader = csv.reader(StringIO.StringIO(content), delimiter=',')
                        linecount = 0
                        for row in reader:
                            if linecount == 0:  # skip the first line of csv file
                                linecount += 1
                                continue
                            linestring = ''.join(row)
                            column = []
                            for col in row:
                                column.append(col)
                            if linestring.startswith('Best value'):
                                break

                            # obtains each iteration result found by a specific execution of a specific node (can be parallel)
                            iteration = long(linecount)
                            value = float(column[1])
                            pos_value = float(column[2])
                            neg_value = float(column[3])
                            K = long(column[4])
                            time = float(column[5])

                            key = file_list[count]
                            iter_data = [iteration, value, pos_value, neg_value, K, time, key]

                            if not full_iteration_data.has_key(key):
                                full_iteration_data[key] = [iter_data]
                            full_iteration_data[key].append(iter_data)
                            linecount += 1
                        count = count - 1
                    finally:
                        content_file.close()

                    # for each execution of the algorithm
                    for executionId in sorted(full_iteration_data.iterkeys()):
                        executionIdIterations = full_iteration_data[executionId]
                        desired_iteration = -1
                        # for each iteration of an execution of the algorithm
                        for iter_data in executionIdIterations:
                            # find out in which iteration the algorithm reaches the target solution value
                            sol_value = float(iter_data[1])
                            if sol_value <= target_value:
                                desired_iteration = iter_data
                                #print "Found target at " + iter_data[5] + " seconds.\n"
                                break
                        # stores the info about the target iteration when the target solution was found
                        if desired_iteration < 0:
                            # if the code gets to this point, it means the algorithm was running in parallel and one
                            # of the nodes did not get to the best solution => can be safely ignored
                            continue
                        jobId = executionId[:executionId.rfind('/')]
                        jobId = jobId[jobId.rfind('/') + 1:]
                        if target_iteration_history.has_key(jobId):
                            current_iteration_info = target_iteration_history[jobId]
                            current_iteration_time = float(current_iteration_info[5])
                            new_iteration_time = float(desired_iteration[5])
                            if current_iteration_time > new_iteration_time:
                                # replaces current iteration info with another one, which is faster
                                # only happens in parallel algorithms
                                target_iteration_history[jobId] = desired_iteration
                                # print "found better result in the same job: " + str(jobId) + ", time = " + str(new_iteration_time) +  "\n"
                        else:
                            target_iteration_history[jobId] = desired_iteration
                # end process all result files for one execution of the algorithm
        # end process all executions of the algorithm

    # now we know the time-to-target of each algorithm execution for this instance
    print "Total: " + str(len(target_iteration_history)) + " results"
    print "Time, Iter, Info"
    tttlist = []
    for key in target_iteration_history.iterkeys():
        iteration_info = target_iteration_history[key]
        tttlist.append(iteration_info[5])
        print "%s %s %s" % (str(iteration_info[5]), str(long(iteration_info[0])), str(iteration_info))

    print "\nFinished TTT processing for instance " + str(filename) + ".\n"
    return tttlist


def generate_box_plot(data_to_plot, instance_names, labels, result_file_prefix):
    # group boxplots - http://stackoverflow.com/questions/20365122/how-to-make-a-grouped-boxplot-graph-in-matplotlib
    fig, axes = plt.subplots(nrows=len(instance_names), sharex=True, figsize=(20, 27)) #ncols=len(labels)) #, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # Create a figure instance
    #fig = plt.figure(1, figsize=(30, 20))

    # Create an axes instance
    #ax = fig.add_subplot(111)

    # Create the boxplot
    ## add patch_artist=True option to ax.boxplot()
    ## to get fill color
    #bp = ax.boxplot(data_to_plot, patch_artist=True, vert=False) #, showfliers=False)

    for ax, name in zip(axes, sorted(instance_names)):
        print "Processing instance " + name

        bp = ax.boxplot([data_to_plot[name][item] for item in xrange(0, len(labels))], patch_artist=True, vert=False) #, showfliers=False)
        ax.set(yticklabels=labels) #, ylabel=name)
        # remove the .g file extension from instance name
        if '.g' in name:
            name = name[:name.rfind('.g')]
        if name.find('file_') >= 0:  # remove file_ prefix from graph files and replace it with 'n='
            name = name[5:]
            n = name[:name.find('_')]
            d = name[name.find('_')+1:name.rfind('_')]
            d_minus = name[name.rfind('_')+1:]
            label_name_1 = 'n = ' + n
            label_name_2 = 'd = ' + d
            label_name_3 = 'd- = ' + d_minus
        else:  # chinese instance files
            c = name[name.find('c')+1:name.find('n')]
            n = name[name.find('n')+1:name.find('k')]
            k = name[name.find('k')+1:name.find('pin')]
            pin = name[name.find('pin')+3:name.find('p-')]
            p_minus = name[name.find('p-')+2:name.find('p+')]
            p_plus = name[name.find('p+')+2:]
            label_name_1 = 'c = ' + c + ', n = ' + n + '            '
            label_name_2 = 'k = ' + k + ', pin = ' + pin + '          '
            label_name_3 = 'p- = ' + p_minus + ', p+ = ' + p_plus + '        '
        ax.set_ylabel(label_name_1 + '\n' + label_name_2 + '\n' + label_name_3 + '\n\n', rotation = 0, labelpad = 70)
        ax.set_xlabel("Execution time (s)", labelpad=20)
        ax.tick_params(axis='y', which='major', pad=5)

        #ax.margins(0.05)  # Optional

        ## change outline color, fill color and linewidth of the boxes
        count = 0
        for box in bp['boxes']:
            # change outline color
            if count % 2 == 0:
                box.set(color='#DAA520', linewidth=2)
            else:
                box.set(color='#7570b3', linewidth=2)
            # change fill color
            if count % 2 == 0:
                #box.set(facecolor='#1b9e77')
                box.set(facecolor='#FFFF66')
            else:
                box.set(facecolor='#800080')
            count += 1

        ## change color and linewidth of the whiskers
        for whisker in bp['whiskers']:
            whisker.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='#7570b3', linewidth=2)

        ## change color and linewidth of the medians
        count = 0
        for median in bp['medians']:
            if count % 2 == 0:
                median.set(color='#654321', linewidth=2)
            else:
                median.set(color='#b2df8a', linewidth=2)
            count += 1

        ## change the style of fliers and their fill
        for flier in bp['fliers']:
            flier.set(marker='o', color='#e7298a', alpha=0.5)

        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

    ## Custom y-axis labels
    #ax.set_yticklabels(instance_names)

    # Save the figure
    plt.tight_layout()
    plt.subplots_adjust(hspace=.01)

    fig.show()
    print "Saving box plot png file to " + result_file_prefix + '-box_plot.png'
    fig.savefig(result_file_prefix + '-box_plot.png') #, bbox_inches='tight')
    pp = PdfPages(result_file_prefix + '-box_plot.pdf')
    pp.savefig(plt.gcf())
    pp.close()


# Originally in http://stackoverflow.com/questions/11882393/matplotlib-disregard-outliers-when-plotting
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


if __name__ == "__main__":
    main(sys.argv[1:])

