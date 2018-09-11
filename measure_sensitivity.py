#!/usr/bin/python

# analyze all reads, and determine how many and which overlaps should be reported
# analyze reported overlaps
# compare ideal and reported overlaps

import re
from operator import itemgetter

# parse and extract integers
def parse(line):
	return map(int, re.findall('\d+', line))


# analyze all reads
# data header: id, startpos, original length, erronous length
f1 = open('reference.fasta','r')
f2 = open('reads.fasta','r')
list1 = []		# contains info from f1
list2 = []		# contains info from f2

for line in f1:
	if ">" in line:
		#print(line)
		t = parse(line)
		list1.append(t)

for line in f2:
	if ">" in line:
		#print(line)
		t = parse(line)
		list2.append(t)

f1.close()
f2.close()

### cal max diff between perfect read length and erronous read length
"""list3 = []
max1 = 0
for t in list1:
	max1 = max(max1, float(abs(t[3]-t[2]))/t[3])

print(max1)
exit()#"""

all_overlaps = []		# contains all theoretical overlaps

## all positions wrt original genome
## a1: startpos of read in ref.fa, a2: endpos of read in ref.fa
## b1: startpos of read in reads.fa, b2: endpos of read in reads.fa
## o1: startpos of overlap, o2: endpos of overlap
## sax: position of start end stop of subread, wrt the read
##		margin_right is larger, because the right is more flexible due to errors
max_length = 0
margin_left = 100		# amount of bases with which the theoretical overlap is extended (if possible)
margin_right = 0.15		# observed max diff between original and erronous read length
for (idx1,r1) in enumerate(list1):
	for (idx2,r2) in enumerate(list2):
		a1 = r1[1]
		b1 = r2[1]
		a2 = a1 + r1[2]
		b2 = b1 + r2[2]
		if a2 < b1 or b2 < a1:
			continue
		o1 = max(a1, b1)
		o2 = min(a2, b2)
		ovl_length = o2 - o1
		all_overlaps.append([idx1,idx2,a1,a2,b1,b2,o1,o2,ovl_length])
		max_length = max(max_length, ovl_length)
		print("r1: %d, r2: %d, a1: %d, a2: %d, b1: %d, b2: %d, o1: %d, o2: %d, length: %d" % (idx1, idx2, a1, a2, b1, b2, o1, o2, ovl_length))
		sa1 = max(a1, o1 - margin_left) - a1
		sa2 = min(min(a2, o2 + int(margin_right*ovl_length)) - a1, list1[idx1][3])
		sb1 = max(b1, o1 - margin_left) - b1
		sb2 = min(min(b2, o2 + int(margin_right*ovl_length)) - b1, list2[idx2][3])
		print("r1: %d, r2: %d, sa1: %d, sa2: %d, sb1: %d, sb2: %d" % (idx1, idx2, sa1, sa2, sb1, sb2))
		print()
print("Max overlap length: " + str(max_length))














# analyze reported overlaps
"""f1 = open('out.darwin')
for line in f1:
	ovl = parse(line)
	print(ovl)

f1.close()
"""
