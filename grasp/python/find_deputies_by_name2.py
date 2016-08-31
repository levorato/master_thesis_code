# script para encontrar os ids dos deputados a partir de uma lista de seus nomes

# to run this script, please install:
# sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose python-tk

import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import re
import math
import pandas as pd

import HTML


def main(argv):

    csv.field_size_limit(1000000000)

    folder = ''
    filter = 'deputados.csv'
    deputies_name_file = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["folder=", "filter=", "depfile="])
    except getopt.GetoptError:
        print 'find_deputies_by_name.py --folder <folder> --filter <filter> --depfile <deputies_name_file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'find_deputies_by_name.py --folder <folder> --filter <filter> --depfile <deputies_name_file>'
            sys.exit()
        elif opt in ("-i", "--folder"):
            folder = arg
        if opt in ("-f", "--filter"):
            filter = arg
        if opt in ("-p", "--depfile"):
            deputies_name_file = arg

    if folder == '':
        print 'Please specify the deputies ids csv dir'
        sys.exit()
    if deputies_name_file == '':
        print 'Please specify the deputies names list path'
        sys.exit()
    print 'Input dir (for csv id files) is ', folder
    print 'File filter is ', filter
    print 'deputies name list file is ', deputies_name_file

    all_deputies_id_df = pd.DataFrame({'A': []})

    for root, subFolders, files in os.walk(folder):
        # sort dirs and files
        subFolders.sort()
        files.sort()
        print "Processing folder " + ''.join(root)
        if len(files):
            file_list = []
            file_list.extend(glob.glob(os.path.join(root, "*" + filter)))
            file_count = len(file_list) - 1

            while file_count >= 0:
                print "Processing csv file " + file_list[file_count] + "\n"

                deputies_id_df = pd.read_csv(file_list[file_count], encoding="utf-8-sig", sep=';', header=0, index_col=0)
                                      # names=['id_dep', 'id_vertex', 'name', 'party', 'state'])

                if all_deputies_id_df.empty:
                    all_deputies_id_df = deputies_id_df
                else:
                    all_deputies_id_df = all_deputies_id_df.append(deputies_id_df)

                file_count = file_count - 1
            # end process results
    # end for
    print "\nSuccessfully processed all deputies id csv files.\n"
    all_deputies_id_df = all_deputies_id_df.reset_index()
    print all_deputies_id_df

    deputies_names_list_df = pd.read_csv(deputies_name_file, encoding="utf-8-sig", sep=';', header=0)
    #deputies_names_list_df['name_upper'] = deputies_names_list_df[str(deputies_names_list_df.columns[0])].str.upper()
    print deputies_names_list_df
    print str(deputies_names_list_df.columns[2])
    deputies_names_list_df['name_upper'] = deputies_names_list_df[str(deputies_names_list_df.columns[2])].str.upper()
    #deputies_names_list_df['name_upper'] = deputies_names_list_df['name_upper'].apply(lambda x: u', '.join(x).encode('utf-8').strip())
    print deputies_names_list_df
    print 'Tamanho da lista de deputados: '
    print len(deputies_names_list_df.index)
    name_and_deputies_id_df = pd.merge(deputies_names_list_df, all_deputies_id_df, left_on = ['name_upper'], right_on=['nomeParlamentar'], how='left')
    print len(name_and_deputies_id_df.index)
    print name_and_deputies_id_df
    name_and_deputies_id_df = name_and_deputies_id_df.drop_duplicates(subset = 'Partido')
    #print name_and_deputies_id_df
    #name_and_deputies_id_df = name_and_deputies_id_df['Partido', 'idCadastro']

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    dest_filename = deputies_name_file[deputies_name_file.rfind(os.path.sep)+1:deputies_name_file.rfind('.')] + '_with_ids.xlsx'
    writer = pd.ExcelWriter(os.path.join(folder, dest_filename), engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    name_and_deputies_id_df.to_excel(writer, sheet_name='ids')

    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
    print "\nExported deputies id list file to Excel file.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
