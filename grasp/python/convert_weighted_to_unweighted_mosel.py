import sys, getopt
import csv
import StringIO
import glob
import os
import os.path
import argparse


def main(argv):
    csv.field_size_limit(1000000000)

    parser = argparse.ArgumentParser(description='Convert mosel graph files (.mos/.g) to unweighted.')
    parser.add_argument('folders', nargs='+',
                        help='the folders containing the mos files')
    parser.add_argument('--filefilter', default='*.g', required=False,
                        help='the file extension for mosel graph files (default: *.g)')
    args = parser.parse_args()
    folders = args.folders
    filter = args.filefilter

    args = parser.parse_args()

    print 'Input folders are ', folders
    print 'File filter is ', filter

    processMoselFiles(folders, filter)


def processMoselFiles(folders, filter):
    for folder in folders:
        print "Processing folder " + ''.join(folder)

        for root, subFolders, files in os.walk(folder):
            # sort dirs and files
            subFolders.sort()
            files.sort()

            if (len(files)):
                file_list = []
                file_list.extend(glob.glob(root + "/" + str(filter)))
                count = len(file_list) - 1
                print 'Found ' + str(count) + ' files'
                directory = root + '/unweighted'
                if not os.path.exists(directory):
                    os.makedirs(directory)

                while count >= 0:
                    # Process all files from all algorithm executions, including those in parallel
                    # print "Processing file " + file_list[count] + "\n"
                    filename = file_list[count]
                    if '/' in filename:  # linux
                        filename = filename[filename.rfind('/') + 1:]
                        linux = True
                    else:  # windows
                        filename = filename[filename.rfind('\\') + 1:]
                        linux = False
                    content_file = open(file_list[count], 'r')
                    if linux:
                        output_file = open(directory + '/' + filename, 'w')
                    else:
                        output_file = open(directory + '\\' + filename, 'w')
                    try:
                        content = content_file.read()

                        reader = csv.reader(StringIO.StringIO(content), delimiter=')')
                        for row in reader:
                            linestring = ''.join(row)
                            column = []
                            #print linestring
                            #print row
                            for col in row:
                                column.append(col)
                            if len(row) == 0:
                                continue

                            if "(" in linestring:
                                weight = float(column[1])
                                if weight > 0:
                                    s_weight = '1'
                                else:
                                    s_weight = '-1'
                                output_file.write(column[0] + ')' + s_weight + '\n')
                            else:
                                output_file.write(column[0] + '\n')
                        count = count - 1
                    finally:
                        content_file.close()
                        output_file.close()

                        # end loop
                        # process last file
        print "\nDone.\n"


if __name__ == "__main__":
    main(sys.argv[1:])
