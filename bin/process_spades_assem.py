#!/usr/bin/env python

# MIT License
#
# Copyright (c) Alexandra Weisberg
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# Parse inputs

#Give contigs file, prints contigs larger than 1k
from Bio import SeqIO
import sys
import argparse
import os

def extant_file(x):
    """
    'Type' for argparse - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        # Argparse uses the ArgumentTypeError to give a rejection message like:
        # error: argument input: x does not exist
        raise argparse.ArgumentTypeError("{0} does not exist".format(x))
    return x

def median(lst):
    lst = sorted(lst)
    n = len(lst)
    if n < 1:
        return None
    if n % 2 == 1:
        return lst[n//2]
    else:
        i = n//2
        return (lst[i - 1] + lst[i]) / 2

parser = argparse.ArgumentParser(description='Extract contigs based on coverage cut-off from FASTA file. Also optionally remove palindromes.')
parser.add_argument('fasta',help='FASTA file with contigs/scaffolds of interest',type=extant_file)
parser.add_argument('--palindrome',help='Cut contigs with palindromic sequence in half',action='store_true',default=False)
parser.add_argument('--summary',help='Filename for a summary file with predicted contig types and palindromic contigs',type=str)
parser.add_argument('--cov_cutoff',help='Keep contigs above this coverage value (default = 5)',type=float,default=5.0)
parser.add_argument('--cov_max',help='Keep contigs below this cov value (default = None)',type=int,default=0)
parser.add_argument('--len_cutoff',help='Keep contigs above this length value (default = 500)',type=int,default=500)
parser.add_argument('--len_max',help='Keep contigs below this length value (default = None)',type=int,default=0)
parser.add_argument('--verbose',help='Print progress messages (default = off)',action='store_true',default=False)

args = parser.parse_args()

summary = {}
coverage = {}
palindromes = {}

infile = SeqIO.parse(args.fasta, 'fasta')

for rec in infile:
    if args.verbose: sys.stderr.write("Working on {}\n".format(rec.id))
    header = rec.id
    data = header.split('_')
    try:
        float(data[5])
    except IndexError:
        sys.stderr.write("Unable to find coverage for sequence {}\n".format(rec.id))
    except ValueError:
        sys.stderr.write("Unable to find valid coverage from header {} - value {}\n".format(rec.id,data[5]))
    else:
        cov = float(data[5])
        if cov > args.cov_cutoff and len(rec.seq) > args.len_cutoff:
            if args.len_max > 0 and len(rec.seq) > args.len_max:
                continue
            else:
                pass
            if args.cov_max > 0 and cov > args.cov_max:
                continue
            else:
                pass
            coverage[rec.id] = int(cov)
            if args.palindrome:
                seq = rec.seq
                revcom = rec.seq.reverse_complement()
                count = 0
                for a, b in zip(seq,revcom):
                    if a != b:
                        count += 1
                    if count > 300:
                        break
                if count < 100:
                    summary[rec.id] = 'palindrome'
                    data[3] = int(len(seq)/2)
                    rec.seq = rec.seq[0:data[3]]
                    data[3] = str(data[3])
                    data.extend(["palindrome","clipped"])
                    rec.id = "_".join(data)
                palindromes[rec.id] = str(count)

            if 'palindrome' in rec.id:
                summary[rec.id] = 'palindrome'
            SeqIO.write(rec,sys.stdout,'fasta')

if args.palindrome:
    try:
        fh = open("palindrome.stats", "w")
    except IOError:
        sys.stderr.write("Unable to open file palindrome.stats for writing. Check path and permissions and try again.\n")
    else:
        for contig in palindromes:
            if palindromes[contig] == 301:
                fh.write("\t".join([contig,"NOT"]) + "\n")
            else:
                fh.write("\t".join([contig,palindromes[contig]]) + "\n")
        fh.close()

if args.summary:
    try:
        fh = open(args.summary, "w")
    except IOError:
        sys.stderr.write("Unable to open file {} for writing. Check path and permissions and try again.\n".format(args.summary))
    else:
        if len(coverage) == 0:
            mid = 0
        else:
            mid = median(coverage.values())
        mid = mid*.66
        for contig in coverage:
            if coverage[contig] < mid:
                summary[contig] = 'plasmid'
        for contig in summary:
            fh.write("\t".join([contig,summary[contig]]) + "\n")
        fh.close()

