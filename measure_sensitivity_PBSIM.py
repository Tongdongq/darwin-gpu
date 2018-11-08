#!/usr/bin/python

# analyze all reads, and determine how many and which overlaps should be reported
# analyze reported overlaps
# compare ideal and reported overlaps

import re, subprocess, bisect, sys, time, os
from operator import itemgetter

# parse and extract integers
def parse(line):
	return map(int, re.findall('\d+', line))

daligner = 0
ref = 0
extra = 1 # if extra is 1, add the reverse of the AB overlap (BA), should increase sensitivity, and halve specificity
remove_trivial = 1 # if equal to 1, removes all 'read X ovls with read X' overlaps
recompute_tovls = 0 # if equal to 1, recomputes all theoretical overlaps based on the input files, otherwise, look for w_theoretical_overlaps

# thresholds to filter heuristic overlaps
score_thres = 800
min_length = 990

if len(sys.argv) > 1:
	arg = sys.argv[1]
	if arg == 'daligner':
		print("Heuristic is Daligner")
		daligner = 1
	if arg == 'ref' or arg == 'reference':
		print("Reference-based")
		ref = 1
if len(sys.argv) > 2:
	arg = sys.argv[2]
	if arg == 'daligner':
		print("Heuristic is Daligner")
		daligner = 1
	if arg == 'ref' or arg == 'reference':
		print("Reference-based")
		ref = 1

if daligner == 0:
	print("Heuristic is Darwin")
if ref == 0:
	print("De-novo based")

# analyze all reads
# data header: id, startpos in genome, length of overlap
if daligner == 0:
	f1 = open('reference.fasta','r')
	f2 = open('reads.fasta','r')
else:
	f1 = open('../PBSIM/src/sd_0001.fasta','r')
	f2 = open('../PBSIM/src/sd_0002.fasta','r')
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

print("Parsed fasta files")
print("Num reads: %d %d" % (len(list1), len(list2)))

max_length = 0		# max length of simulated read

if ref == 0:
	for read in list1:
		max_length = max(max_length, read[2])
for read in list2:
	max_length = max(max_length, read[2])

print("Max length simulated read: %d" % max_length)

all_theoretical_overlaps = []

## all positions wrt original genome
## a1: startpos of read in ref.fa, a2: endpos of read in ref.fa
## b1: startpos of read in reads.fa, b2: endpos of read in reads.fa
## o1: startpos of overlap, o2: endpos of overlap
## sax: position of start end stop of subread, wrt the read
## 0: is added to measure sensitivity

## another, slower way to find theoretical ovls:
### - sort both lists on startpos in genome
### - for each read in list1:
### -- start at the first read that has start_pos2 > start_pos1 - max_length
### -- check for ovls
### -- stop at the first read that has start_pos2 > start_pos1 + max_length
## this approach took 2m17 to find tovls, as oppose to 1m29

if ref == 0:
	if recompute_tovls == 1:
		for (idx1, r1) in enumerate(list1):
			if idx1 % 4000 == 0:
				print('idx1: %d' % idx1)
			for (idx2, r2) in enumerate(list2):
				a1 = r1[1]
				b1 = r2[1]
				a2 = a1 + r1[2]
				b2 = b1 + r2[2]
				if a2 < b1 or b2 < a1:
					continue
				o1 = max(a1, b1)
				o2 = min(a2, b2)
				ovl_length = o2 - o1
				#print("r1: %d, r2: %d, a1: %d, a2: %d, b1: %d, b2: %d, o1: %d, o2: %d, length: %d" % (idx1, idx2, a1, a2, b1, b2, o1, o2, ovl_length))
				sa1 = o1 - a1
				sa2 = o2 - a1
				sb1 = o1 - b1
				sb2 = o2 - b1
				#print("%d %d" % (sa2-sa1, sb2-sb1))
				if ovl_length >= 1000:
					all_theoretical_overlaps.append((idx1, idx2, sa1, sa2, sb1, sb2, 0))
				#print("r1: %d, r2: %d, sa1: %d, sa2: %d, sb1: %d, sb2: %d" % (idx1, idx2, sa1, sa2, sb1, sb2))
				#print()

		print("Num theoretical ovls: %d" % len(all_theoretical_overlaps))

		fout = open('w_theoretical_ovls','w')

		for tovl in all_theoretical_overlaps:
			fout.write(str(tovl))
			fout.write('\n')

		fout.close()
	# above part to calculate tovls, below part to read precalculated tovls from file
	else:
		print("File with tovls created %s" % time.ctime(os.path.getmtime('w_theoretical_ovls')))
		fout = open('w_theoretical_ovls','r')
		for line in fout:
			all_theoretical_overlaps.append(parse(line))
		fout.close()
		print("Read tovls from file")
		print("Num theoretical ovls: %d" % len(all_theoretical_overlaps))

if remove_trivial == 1:
	all_theoretical_overlaps = [ovl for ovl in all_theoretical_overlaps if ovl[0] != ovl[1]]
	print("Num non-trivial theoretical ovls: %d" % len(all_theoretical_overlaps))

# analyze reported overlaps by heuristic aligner
## darwin:
## ref_id, gen_pos, ovl_length, \
## read_id, gen_pos, ovl_length, \
## ref_start, ref_end, read_start, read_end, score, comp, 0
## daligner:
## ref_id, read_id, \
## ref_start, ref_end, \
## read_start, read_end, 0
## last 0 is used to find False Positives (FP)
all_heuristic_overlaps = []
if daligner == 0:
	print('Reading darwin overlaps')
	f1 = open('out.darwin')
	for line in f1:
		l = parse(line)
		if len(l) < 10:
			print('WARNING this darwin overlap does not have enough information')
			print(line)
			print(l)
		l.append(0)
		all_heuristic_overlaps.append(l)
		if extra == 1:
			l = [l[3],l[4],l[5],l[0],l[1],l[2],l[8],l[9],l[6],l[7],l[10],l[11],0]
			all_heuristic_overlaps.append(l)
	f1.close()
else:
	print('Reading daligner overlaps')
	# daligner must have run already, with './run.sh 8 | tee >(grep "ids" > tmp) | grep -v "ids"'
	# only consider ref_id and read_id for now, for coordinates, the commas need to be removed
	f1 = open('../original/DALIGNER-old/tmp')
	for line in f1:
		l = parse(line)
		if len(l) != 7:
			print('WARNING this daligner overlap does not have enough information')
			print(line)
			print(l)
		l.append(0)
		all_heuristic_overlaps.append(l)
	f1.close()

print("Num heuristic overlaps: %d" % len(all_heuristic_overlaps))
if remove_trivial == 1:
	all_heuristic_overlaps = [ovl for ovl in all_heuristic_overlaps if ovl[0] != ovl[1]]
	print("Num non-trivial heuristic overlaps: %d" % len(all_heuristic_overlaps))

# filter out some heuristic overlaps
## example criteria: length of overlap, score

hidx = 3			# idx where the read_id is, inside the heuristic ovl
last_idx = 12		# idx where the extra 0 is placed, which is used to count FP
score_idx = 10		# idx where the score is, inside the heuristic ovl
sa1 = 6; sa2 = 7; sb1 = 8; sb2 = 9

if ref == 1:
	min_length = 0
	if daligner == 1:
		last_idx = 999
		sa1 = 2; sa2 = 3; sb1 = 4; sb2 = 5; score_idx = -1
	else:
		last_idx = 10
		sa1 = 4; sa2 = 5; sb1 = 6; sb2 = 7; score_idx = 8
else:
	if daligner == 1:
		hidx = 1
		last_idx = 7
		sa1 = 2; sa2 = 3; sb1 = 4; sb2 = 5; score_idx = 6


all_heuristic_overlaps = [ovl for ovl in all_heuristic_overlaps if ovl[sa2]-ovl[sa1] >= min_length and ovl[sb2]-ovl[sb1] >= min_length]
if score_idx == -1:
	# daligner score = length ovl - diffs
	all_heuristic_overlaps = [ovl for ovl in all_heuristic_overlaps if ovl[sa2] - ovl[sa1] - ovl[6] >= score_thres]
else:
	all_heuristic_overlaps = [ovl for ovl in all_heuristic_overlaps if ovl[score_idx] >= score_thres]

print("Filter: score_thres: %d, min_length: %d" % (score_thres, min_length))
print("Num heuristic overlaps after filtering: %d" % len(all_heuristic_overlaps))

if len(all_heuristic_overlaps) == 0:
	print("WARNING no heuristic overlaps left, filter is probably too strict")
	exit()

# compare theoretical overlaps and heuristic overlaps

## theoretical ovl: [idx1, idx2, sa1, sa2, sb1, sb2, 0]
## compare theoretical and heuristic overlaps
FN = 0				# false negatives, a tovl has no matching hovl, thus heuristic missed one
FP = 0				# false positives, a hovl has no matching tovl, thus should not exist
TP = 0				# true positives
all_heuristic_overlaps.sort(key=itemgetter(0))

# get ordered list of ref_id from hovls
tmp_list = zip(*all_heuristic_overlaps)[0]

if ref == 0:
	for tovl in all_theoretical_overlaps:
		fn = 1
		# start searching for matches where the ref_id matches
		idx = bisect.bisect_left(tmp_list, tovl[0])
		#print("ids %d %d, start idx: %d" % (tovl[0], tovl[1], idx))
		while idx < len(all_heuristic_overlaps) and all_heuristic_overlaps[idx][0] == tovl[0]:
			hovl = all_heuristic_overlaps[idx]
			if tovl[0] == hovl[0] and tovl[1] == hovl[hidx]:
				fn = 0
				TP = TP + 1
				hovl[last_idx] = 1
				#print tovl
				# remove some info from darwin overlap
				#movl = [y for x in [[hovl[0]], [hovl[3]], hovl[6:]] for y in x]
				#print movl
				#print()

				#if n % 10 == 0:
					#print("ref_id, read_id, ref_start, ref_end, read_start, read_end, score")
			idx += 1
		if fn == 1:
			FN = FN + 1
			#if FN < 500 and FN > 490:
			#	print tovl
	for hovl in all_heuristic_overlaps:
		if hovl[last_idx] == 0:
			FP = FP + 1
			#print hovl
else:
	# fp is a list
	## 0: read not mapped
	## 1: read is mapped to the correct spot
	fp = [0] * len(list2)
	unique_hovls = []								# contains the hovl with the highest score, for each readID
	all_heuristic_overlaps.sort(key=itemgetter(1))	# sort by readID
	readID = all_heuristic_overlaps[0][1]
	maxScore = -1
	maxScoreIdx = 0
	for i,h in enumerate(all_heuristic_overlaps):
		if h[1] == readID:
			if daligner == 1:
				score = h[sa2]-h[sa1]-h[6]
			else:
				score = h[8]
			if score > maxScore:
				maxScore = score
				maxScoreIdx = i
		else:
			t = all_heuristic_overlaps[maxScoreIdx]
			#print("ID: %d, idx: %d, %s, %s" % (readID, maxScoreIdx, str(t), str(h)))
			unique_hovls.append(all_heuristic_overlaps[maxScoreIdx])
			readID = h[1]
			maxScoreIdx = i
			maxScore = -1

	print("Hovls with unique readID: %d" % len(unique_hovls))#"""

	for hovl in unique_hovls:
	#for hovl in all_heuristic_overlaps:
		if daligner == 0:
			read_id = hovl[1]
			gen_pos = hovl[2]
			ref_start = hovl[4]
			ref_end = hovl[5]
		else:
			if len(hovl) != 8:
				print("ERROR daligner too little info")
				print hovl
			read_id = hovl[1]
			gen_pos = list2[read_id][1]
			ref_start = hovl[2]
			ref_end = hovl[3]
		if gen_pos > ref_start-50 and gen_pos < ref_start+50:
			TP += 1
		else:
			FP += 1
			#if FP < 50:
			#	print hovl
		fp[read_id] = 1

	FN = len(list2) - sum(fp)

	"""i = 0
	for (idx,f) in enumerate(fp):
		if f == 0:
			if i > 100:
				print list2[idx]
			i += 1
		if i > 110:
			break#"""


print("TP: %d" % TP)
print("FN: %d" % FN)
print("FP: %d" % FP)

print("sensitivity: %f" % (float(TP)/(TP+FN)))
print("specificity: %f" % (float(TP)/(TP+FP)))


