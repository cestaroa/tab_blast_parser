#!/usr/bin/python
import os.path
import sys
import re
import argparse
#
#tabular blast, give back a dictionary with query_id as key and an array of dicts
#
def build_by_field(my_file,tab_field) :
    by_query={}
    blast_file=open(my_file)
    my_test=False
    my_test_len=False
    #also filter by length ratio
    for l in blast_file.readlines() :
        l=l.rstrip()
        #outfmt 7
        if l[0]=='#' :
            #avoids overwriting fields
            if my_test :
                continue
            #catch field names
            if l[0:8]=='# Fields' :
                l=l.split(', ')
                l[0]=l[0][10:]#remove Fields
                l[2]=l[2][2:]#remove perc sign
                tab_field=l
                my_test=True
        else :
            l=l.split("\t")
            my_dict={}
            for i in range(0,len(tab_field)) :
                my_dict[tab_field[i]]=l[i]
            #
            #fix value types
            for i in tab_field[2:] :
                my_dict[i]=float(my_dict[i])
            #
            #used to filter results
            my_dict['filter']=False
            if my_dict['query id'] in by_query :
                by_query[my_dict['query id']].append(my_dict)
            else :
                by_query[my_dict['query id']]=[my_dict]
    return(by_query,tab_field)
#
def filter_by_field() :
    for hsp in blast_by_query.values() :
        for h in hsp :
            for f in filter_field :
                #print f,h[f],filter_field[f]
                if f=='evalue' :
                    if h[f]>filter_field[f] :
                        h['filter']=True
                else :
                    if h[f]<filter_field[f] :
                        h['filter']=True
#
def get_ratio() :
    for qid in blast_by_query :
        for h in blast_by_query[qid] :
            f1=h['alignment length']/h['query length']
            f2=h['alignment length']/h['subject length']
            f3=h['query length']/h['subject length']
            h['alignment/query']=round(f1,3)
            h['alignment/subject']=round(f2,3)
            h['query/subject']=round(f3,3)
#
def sort_by(my_field_name,my_hsp_ary) :
    tmp=[]
    to_return=[]
    #add to tmp all values found in the array of hsp
    for h in my_hsp_ary :
        tmp.append(h[my_field_name])
    #
    tmp=sorted(set(tmp))
    #for lengths and scores is better to have a reverse sorting
    if my_field_name=='alignment length' or my_field_name=='bit score' :
        tmp=reversed(tmp)
    #
    for t in tmp :
        to_del=[]
        for i in range(0,len(my_hsp_ary)) :
            if my_hsp_ary[i][my_field_name]==t :
                to_return.append(my_hsp_ary[i])
                #to remove already added hsp
                to_del.append(i)
        for i in reversed(sorted(to_del)) :
            del(my_hsp_ary[i])
    return(to_return)
#
def print_by_query(my_sort_par) :
    to_print=[]
    for hsp in blast_by_query.values() :
        for h in sort_by(my_sort_par,hsp) :
            if h['filter'] :
                continue
            tmp=[]
            for t in tab_field :
                tmp.append(str(h[t]))
            to_print.append("\t".join(tmp))
    print "#"+"\t".join(tab_field)
    print "\n".join(to_print)
#
#
#
#general params
#default fields of tabular format
tab_field=['query id','subject id','identity','alignment length','mismatches','gaps','q. start','q. end','s. start','s. end','evalue','bit score']
#
#
#input params
parser = argparse.ArgumentParser(description='Simple parser for blast results in tab format.')
parser.add_argument("blast", type=str, help="tabular blast file")
parser.add_argument("-i", "--identity", type=float, help="identity cut-off value default 30.0",default=30.0)
parser.add_argument("-l", "--length", type=int, help="alignment length cut-off value, default 30",default=30)
parser.add_argument("-e", "--evalue", type=float, help="e-value cut-off, default 0.001",default=0.01)
parser.add_argument("-s","--sort_by", type=str, help="return results sorted by query id and by the indicated feature",default='bit score')
parser.add_argument("-r", "--length_ratios", type=bool, help="add to results the ration among alignemnt len and query/sbjct len (required modified outfmt 7), default False")
args=parser.parse_args()
#
#Parse input params
#
#puts the values for filtering in the dedicated array
filter_field={}
if args.evalue :
    filter_field['evalue']=args.evalue
if args.identity :
    filter_field['identity']=args.identity
if args.length :
    filter_field['alignment length']=args.length
if args.blast :
    if not os.path.isfile(args.blast) :
        print "help"
        sys.exit()
else :
    print "help"
    sys.exit()
#add sorting parameter. Print a message on stderr if parameter does not match with column names
if not args.sort_by in tab_field :
    sys.stderr.write("Warning !!! "+args.sort_by+" is not among the column names. Bit score will be used instead\n")
    args.sort_by='bit score'
else :
    args.sort_by=args.sort_by.strip("'")
#
#main part
#
(blast_by_query,tab_field)=build_by_field(args.blast,tab_field)
#since tab_field array is overwrited by the precedent call, ratio labels has to be added.
if args.length_ratios :
    tab_field+=['alignment/query','alignment/subject','query/subject']
    get_ratio()
#
filter_by_field()
print_by_query(args.sort_by)
