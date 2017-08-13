# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

# <editor-fold desc="Command-line arguments to calculate chinese box plot graph.">
# /home/mlevorato/Downloads/mestrado-output-temp/output/ils-seq-SetParA/agents-seq-l1-ils10s-a0.4-fi-imb-2h
# /home/mlevorato/Downloads/mestrado-output-temp/output/grasp-bullmpi/agents-seq-l1-grasp400s-a1.0-fi-imb-2h
# --labels
# SeqILS
# SeqGRASP
# --excludefiles
# c4n4096k16pin0.8p-0.8p+0.6.g

# --excludefiles
# ml-10M100K_mw0.8000.g-cosine.g
# ml-10M100K_mw0.8000.g-pearson.g
# ml-10M100K_mw0.9000.g-cosine.g
# ml-10M100K_mw0.9000.g-pearson.g
# ml-20m_mw0.6000.g-cosine.g
# ml-20m_mw0.6500.g-cosine.g
# ml-20m_mw0.9000.g-pearson.g
#
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

    parser = argparse.ArgumentParser(description='Compare GRASP / ILS - CC result files, generating box plot graphs based on imbalance value.')
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
    imbalance_file_summary = dict()
    algorithm_list = []
    for i, folder in enumerate(folders):
        print "Processing folder L1 " + ''.join(folder)
        # CC results
        previous_filename = ""
        algorithm_list.append(folder[folder.rfind(os.path.sep) + 1:])

        for root1, subFolders1, files1 in os.walk(folder):
            # sort dirs and files
            subFolders1.sort()
            files1.sort()

            print "Processing folder L2 " + ''.join(root1)
            if (''.join(root1) != folder):  # process only subfolders
                current_folder = (root1[root1.rfind(os.path.sep) + 1 :])
                if '.g' in current_folder:
                    filename = current_folder
                    algorithm = root1[:root1.rfind(filename) - 1]
                    algorithm = algorithm[algorithm.rfind(os.path.sep) + 1:]

                    for root, subFolders, files in os.walk(root1):
                        # Process CC results
                        if os.path.isfile(root + os.path.sep + "cc-result.txt"):
                            filename = (root[:root.rfind(os.path.sep)])
                            filename = filename[filename.rfind(os.path.sep)+1:]
                            datetime = root[root.rfind(os.path.sep) + 1:]
                            skip = False
                            if exclude != None:
                                if filename in exclude:
                                    print "Skipping file " + str(filename)
                                    skip = True
                                for pattern in exclude:
                                    if pattern in filename:
                                        print "Skipping file " + str(filename) + ' due to exclude file pattern'
                                        skip = True
                            if skip:
                                continue
                            best_value = sys.float_info.max

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

                            # captura o resultado de imbalance dadas todas as execucoes de um mesmo grafo / instancia
                            print "Storing result for filename = " + filename + " and algorithm: " + algorithm
                            if not imbalance_file_summary.has_key(filename):
                                imbalance_file_summary[filename] = dict()
                            if not imbalance_file_summary[filename].has_key(algorithm):
                                imbalance_file_summary[filename][algorithm] = []
                            imbalance_file_summary[filename][algorithm].append(float(best_value))

                    # end process all result files of an instance / filename
                # end loop
                # process last file

    data_to_plot = {}
    instance_names = []
    for instance_name in imbalance_file_summary.iterkeys():
        print '\nInstance: ' + str(instance_name) + ':\n'
        result_for_filename = imbalance_file_summary[instance_name]
        instance_names.append(instance_name)
        if not data_to_plot.has_key(instance_name):
            data_to_plot[instance_name] = []
        for algorithm in algorithm_list: #result_for_filename.iterkeys():
            imblist = result_for_filename[algorithm]
            print 'Algorithm: ' + str(algorithm) + ' - total results: ' + str(len(imblist))
            # uses the imbalance series to generate the box plots
            ## combine these different collections into a list
            # Remove the outliers from the list
            imbarray = np.array(imblist)
            filtered = imbarray[~is_outlier(imbarray)]
            print 'Total items after outliar filtering: ' + str(len(filtered)) + '\n'
            data_to_plot[instance_name].append(filtered)
        # end for
    # end for
    experiment_name = ''
    for folder in folders:
        experiment_name = experiment_name + folder[folder.rfind('/') + 1:] + '-vs-'
    tdir = './compare_results'
    if not os.path.exists(tdir):
        os.makedirs(tdir)
    result_file_prefix = tdir + '/' + str(experiment_name)

    print "Generating box plot graph for imbalance series for each instance..."
    generate_box_plot_horizontal(data_to_plot, instance_names, labels, result_file_prefix)
    #generate_box_plot(data_to_plot, instance_names, labels, result_file_prefix)
    print "Box plots successfully generated."



def generate_box_plot_horizontal(data_to_plot, instance_names, labels, result_file_prefix):
    # Add a dummy subplot as first subplot to help correcting ylabel spacing issues (1)
    dummy_instance = "dummy"
    instance_names = sorted(instance_names) + [dummy_instance]
    data_to_plot[dummy_instance] = []
    data_to_plot[dummy_instance].append([])
    data_to_plot[dummy_instance].append([])

    # group boxplots - http://stackoverflow.com/questions/20365122/how-to-make-a-grouped-boxplot-graph-in-matplotlib
    height = 2 * len(instance_names)
    padding = 40  #60
    bbox_to_anchor = (0.98, 0.965)
    if instance_names[0].rfind('file_') >= 0:
        height = 6 #20
        padding = 40
        bbox_to_anchor = (0.98, 0.955)
    elif instance_names[0].find('ml-') >= 0:
        height = len(instance_names)
    elif instance_names[0][0] == 'c':
        height = 1.2 * len(instance_names)
    fig, axes = plt.subplots(nrows=len(instance_names), sharex=False,
                             figsize=(12, height))  # ncols=len(labels)) #, sharex=True, sharey=True)
    if len(instance_names) == 1:
        axes = [axes]

    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # Create a figure instance
    #fig = plt.figure(1, figsize=(30, 20))

    # draw temporary colored lines and use them to create a legend
    hB, = plt.plot([1, 1], '#800080')
    hR, = plt.plot([1, 1], '#DAA520')
    if len(labels) == 1:
        fig.legend((hB, hR), ('(1) ' + labels[0], '(2) ' + labels[0]), loc='upper right', shadow=True,
               bbox_to_anchor=bbox_to_anchor, ncol=2, borderaxespad=0.9)
    else:
        fig.legend((hB, hR), ('(1) ' + labels[0], '(2) ' + labels[1]), loc='upper right', shadow=True,
               bbox_to_anchor=bbox_to_anchor, ncol=2, borderaxespad=0.9)

    # bbox_to_anchor=(0., 1.02, 1., .102))  # bbox_to_anchor=(0, 1))
    hB.set_visible(False)
    hR.set_visible(False)
    # Create an axes instance
    #ax = fig.add_subplot(111)

    ylabels = instance_names
    print ylabels
    axis_count = -1
    for ax, name in zip(axes, instance_names):
        print "Processing instance " + name
        print "Number of series to plot: " + str(len(data_to_plot[name]))
        axis_count += 1

        if len(labels) != len(data_to_plot[name]):
            print 'WARNING: Skipping instance ' + name  + ' due to lack of results!'
            continue
        if len(labels) > 1:
            bp = ax.boxplot([data_to_plot[name][1-item] for item in xrange(0, len(labels))], patch_artist=True, vert=False) #, showfliers=False)
        else:
            bp = ax.boxplot([data_to_plot[name][0]], patch_artist=True, vert=False)  # , showfliers=False)
        if name != 'dummy':
            if len(labels) > 1:
                xmin = min([min(data_to_plot[name][1-item]) for item in xrange(0, len(labels))])
                xmax = max([max(data_to_plot[name][1-item]) for item in xrange(0, len(labels))])
            else:
                xmin = min(data_to_plot[name][0])
                xmax = max(data_to_plot[name][0])
        ax.set_xlim([xmin - 0.01 * xmin, xmax + 0.01 * xmax])
        #axes.set_ylim([ymin, ymax])
        #ax.set(yticklabels=labels) #, ylabel=name)
        ax.set(yticklabels=['(2)', '(1)'])  # , ylabel=name)
        # remove the .g file extension from instance name
        if len(labels) == 1:
            name = ylabels[axis_count]
        else:
            name = ylabels[axis_count - 1]
        additional_line = ''
        if '.g' in name:
            name = name[:name.rfind('.g')]
        if name.find('dummy') >= 0:
            label_name_1 = '\nInstance    '
            label_name_2 = ''
            label_name_3 = '\n'
        elif name.find('file_') >= 0:  # remove file_ prefix from graph files and replace it with 'n='
            name = name[5:]
            n = name[:name.find('_')]
            d = name[name.find('_')+1:name.rfind('_')]
            d_minus = name[name.rfind('_')+1:]
            label_name_1 = 'n = ' + n
            label_name_2 = 'd = ' + d
            label_name_3 = 'd- = ' + d_minus
            additional_line = '\n'
        elif name[0] == 'c':  # chinese instance files
            print 'Chinese random instance type.'
            c = name[name.find('c')+1:name.find('n')]
            n = name[name.find('n')+1:name.find('k')]
            k = name[name.find('k')+1:name.find('pin')]
            pin = name[name.find('pin')+3:name.find('p-')]
            p_minus = name[name.find('p-')+2:name.find('p+')]
            p_plus = name[name.find('p+')+2:]
            label_name_1 = 'c = ' + c + ', n = ' + n + '            '
            label_name_2 = 'k = ' + k + ', pin = ' + pin + '          '
            label_name_3 = 'p- = ' + p_minus + ', p+ = ' + p_plus + '        '
        elif name.find('ml-') >= 0:  # movielens file
            prefix = name[:name.find('_')]
            mw = float(name[name.find('_mw') + 3:name.find('.g')])
            label_name_1 = prefix
            label_name_2 = 'mw = ' + "{0:.2f}".format(mw)
            label_name_3 = ''
        else:
            label_name_1 = '\n' + name
            label_name_2 = ''
            label_name_3 = ''

        if name[0] == 'c':
            ax.set_ylabel('\n' + label_name_1 + '\n' + label_name_2 + '\n' + label_name_3 + '\n\n\n' + additional_line, rotation = 0,  # acrescentar mais um '\n' no final para o texto subir
                      labelpad = padding)
        else:
            ax.set_ylabel('\n\n\n' + label_name_1 + '\n' + label_name_2 + '\n' + label_name_3 + '\n\n' + additional_line, rotation=0, labelpad=padding)
        if axis_count == 0:
            ax.set_xlabel("Solution value                      \n", labelpad=20)
            ax.xaxis.set_label_position('top')
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
        #ax.get_xaxis().tick_bottom()
        #if axis_count == 0:
        ax.get_xaxis().tick_top()
        #if axis_count == 1:


        ax.get_yaxis().tick_left()

    ## Custom y-axis labels
    #ax.set_yticklabels(instance_names)

    # Save the figure
    plt.tight_layout()
    # ajuste da margem entre os subplots
    plt.subplots_adjust(hspace=.8)

    fig.show()
    print "Saving box plot png file to " + result_file_prefix + '-box_plot.png'
    fig.savefig(result_file_prefix + '-box_plot.png') #, bbox_inches='tight')
    pp = PdfPages(result_file_prefix + '-box_plot.pdf')
    pp.savefig(plt.gcf())
    pp.close()


# Vertical box plots
def generate_box_plot(data_to_plot, instance_names, labels, result_file_prefix):

    instance_names = sorted(instance_names)
    # group boxplots - http://stackoverflow.com/questions/20365122/how-to-make-a-grouped-boxplot-graph-in-matplotlib
    width = 12
    height = 6
    padding = 5
    bbox_to_anchor = (0.98, 0.965)
    rotation = 90
    if instance_names[0].rfind('file_') >= 0:
        width = 9
        height = 6 #20
        padding = 5
        bbox_to_anchor = (0.98, 0.955)
        rotation = 60
    fig, axes = plt.subplots(ncols=len(instance_names), sharey=True, figsize=(width, height)) #ncols=len(labels)) #, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0)
    fig.subplots_adjust(hspace=0)

    # Create a figure instance
    #fig = plt.figure(1, figsize=(30, 20))

    # draw temporary colored lines and use them to create a legend
    hB, = plt.plot([1, 1], '#800080')
    hR, = plt.plot([1, 1], '#DAA520')
    fig.legend((hR, hB), ('(1) ' + labels[0], '(2) ' + labels[1]), loc='upper right', shadow=True,
               bbox_to_anchor=bbox_to_anchor, ncol=2, borderaxespad=0.)
    # bbox_to_anchor=(0., 1.02, 1., .102))  # bbox_to_anchor=(0, 1))
    hB.set_visible(False)
    hR.set_visible(False)
    # Create an axes instance
    #ax = fig.add_subplot(111)

    ylabels = instance_names
    print ylabels
    axis_count = -1
    print axes
    print instance_names
    if len(instance_names) == 1:
        axes = [axes]
    for ax, name in zip(axes, instance_names):
        print "Processing instance " + name
        axis_count += 1

        bp = ax.boxplot([data_to_plot[name][item] for item in xrange(0, len(labels))], patch_artist=True, vert=True) #, showfliers=False)
        #ax.set(yticklabels=labels) #, ylabel=name)
        ax.set(xticklabels=['(1)', '(2)'])  # , ylabel=name)
        # remove the .g file extension from instance name
        name = ylabels[axis_count]
        additional_line = ''
        if '.g' in name:
            name = name[:name.rfind('.g')]
        if name.find('dummy') >= 0:
            label_name = 'Instance'
        elif name.find('file_') >= 0:  # remove file_ prefix from graph files and replace it with 'n='
            name = name[5:]
            n = name[:name.find('_')]
            d = name[name.find('_')+1:name.rfind('_')]
            d_minus = name[name.rfind('_')+1:]
            label_name = n + '/' + d + '/' + d_minus
        else:  # chinese instance files
            c = name[name.find('c')+1:name.find('n')]
            n = name[name.find('n')+1:name.find('k')]
            k = name[name.find('k')+1:name.find('pin')]
            pin = name[name.find('pin')+3:name.find('p-')]
            p_minus = name[name.find('p-')+2:name.find('p+')]
            p_plus = name[name.find('p+')+2:]
            label_name = c + '/' + n + '/' + k + '/' + pin + '/' + p_minus + '/' + p_plus
        ax.set_xlabel(label_name, rotation = rotation, labelpad = padding)
        if axis_count == 0:
            ax.set_ylabel("Execution time (s)", labelpad=5)
            #ax.yaxis.set_label_position('top')
        ax.tick_params(axis='y', which='major', pad=1)

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
        #ax.get_xaxis().tick_bottom()
        #if axis_count == 0:
        #    ax.get_xaxis().tick_top()
        #if axis_count == 1:

        ax.get_xaxis().tick_bottom()

    ## Custom y-axis labels
    #ax.set_yticklabels(instance_names)

    # Save the figure
    #plt.tight_layout()
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

