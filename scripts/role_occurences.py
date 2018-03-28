#!/usr/bin/env python2

from __future__ import print_function
import sys, time
stime = time.time()
import argparse
import subprocess

# Instantiate the parser
parser = argparse.ArgumentParser(description='Produces 3-column table of roles, genomes, and features.')

parser.add_argument('-r', type=str,
                    help='File of roles')

parser.add_argument('-g', type=str,
                    help='File of genome IDs')

args = parser.parse_args()

roles_file = args.r
genomes_file = args.g
roles_list = []

with open(roles_file) as roles:
    for line in roles:
        roles_list.append(line.rstrip())


# svc_all_features -i ~/Tmp/foo peg | svc_function_of > ~/Tmp/bar

p1 = subprocess.Popen(["svc_all_features", "-i", genomes_file, "peg"], stdout=subprocess.PIPE)
p2 = subprocess.Popen(["svc_function_of"], stdin=p1.stdout, stdout=subprocess.PIPE)
p3 = subprocess.Popen(["svc_functions_to_roles"], stdin=p2.stdout, stdout=subprocess.PIPE)
p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
output = p3.communicate()[0].splitlines()

for x in output:
    fields = x.split('\t')

    if fields[3] in roles_list:
        print(fields[3], fields[0], fields[1], sep='\t')


