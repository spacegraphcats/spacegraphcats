#!/usr/bin/env python

'''****************************************************************************

*   Program - Convert To GFA format
*   Author - Mayank Pahadia


 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.

****************************************************************************'''

'''****************************************************************************

*   To run the program, pass three arguments with the python script on command line.
*   For example - python convertToGFA.py inputFileName outputFileName kmerSize

*   Logic - It reads through the fasta file with all the unitigs information 
*   and link information and outputs it in the GFA format. 


****************************************************************************'''


import sys
import argparse

def write_segment(name,segment,optional,g,links):
    add = ""
    add += "S\t" #for segment
    add += name #id of segment
    add += "\t"
    add += segment #segment itself
    add += "\t"
    for i in optional: #optional tags
        add+=i
        add+="\t"
    #adding Segment to the file
    g.write(add.strip()+"\n")
    for j in links: #adding all the links of the current segment to the GFA file
        g.write(j)


def main():
    parser = argparse.ArgumentParser(description="Convert a bcalm-generated FASTA to a GFA.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('inputFilename', help='Input FASTA file')
    parser.add_argument('outputFilename', help='Output GFA file')
    parser.add_argument('kmerSize', type=int, help='k-mer length')
    parser.add_argument('-s', '--single-directed', action='store_true',
                        help='Avoid outputting the whole skew-simmetric graph and output only one edge between two nodes',
                        dest='single_directed')

    args = parser.parse_args()

    with open(args.inputFilename) as f:
        #name stores the id of the unitig
        #optional is a list which stores all the optional tags of a segment
        #links stores all the link information about a segment
        name = ""
        optional=[]
        links=[]
        g = open(args.outputFilename,'w')
        #adding Header to the file
        k = int(args.kmerSize)
        g.write('H\tVN:Z:1.0\tks:i:%d\n' %k) # includes the k-mer size
        print("GFA file open")
        #firstLine is for implemetation purpose so that we don't add some garbage value to the output file.
        firstLine = 0
        #segment stores the segment till present, in a fasta file, segment can be on many lines, hence we need to get the whole segment from all the lines
        segment = ""
        for line in f:
            line = line.replace("\n","")
            if(line[0]!=">"):
                #segment might be in more than one line, hence we get the whole segment first, and then put it in the GFA file.
                segment += line
            if(line[0]==">"):
                if(firstLine!=0):#if it's not the firstline in the input file, we store the input in GFA format in the output file
                    write_segment(name,segment,optional,g,links)
                    segment = ""

                firstLine = 1
                #once the previous segment and it's information has been stored, we start the next segment and it's information
                a = line.split(" ")
                name=a[0][1:] #get the id
                optional=[]
                links = []
                #we skip the first value because the first value is ">ID"
                for i in range(1,len(a)):
                    #we need this because the line can end with a space, hence we get one extra value in our list.
                    if(a[i]==""):
                        continue
                    if(a[i][0:2] == "MA"): #previous bcalm2 versions had "MA=[xxx]" optional tag as well, kept it just for compatibility, and reformated
                        optional.append(a[i][0:2]+":f:"+a[i][2:])
                    elif(a[i][0:2] == "L:"): #for links
                        b = a[i].split(":")
                        k1 = int(args.kmerSize)-1
                        if args.single_directed:
                            if name < b[2]:
                                links.append("L\t"+name+"\t"+b[1]+"\t"+b[2]+"\t"+b[3]+"\t"+str(k1)+"M\n")
                            elif name == b[2] and not (b[1] == b[3] == '-'): # manage links between the same unitig
                                links.append("L\t"+name+"\t"+b[1]+"\t"+b[2]+"\t"+b[3]+"\t"+str(k1)+"M\n")
                        else:
                            links.append("L\t"+name+"\t"+b[1]+"\t"+b[2]+"\t"+b[3]+"\t"+str(k1)+"M\n")
                    else: #all the other optional tags
                        optional.append(a[i])


        #we will miss the last one, because it won't go into the if condition - if(line[0]==">") and hence won't add the segment to the file.
        write_segment(name,segment,optional,g,links)
        print("done")
        g.close()

if __name__ == "__main__":
    main()
