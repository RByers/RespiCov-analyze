import os
import numpy as np
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from Bio.Data import IUPACData
from IPython.display import clear_output 
from dataclasses import dataclass
import random
import gzip
import time
import json


# Look for at least an 80% match against primers
# Below about 70% we seem to get huge numbers of matches just by chance
# The Random primer gets a 76% match.
MATCH_THRESHOLD = 0.80

# Hits which overlap atleast 80% of a higher-score hit aren't reported
OVERLAP_THRESHOLD = 0.80

# Use a fixed random seed for reproducibility
rnd = random.Random(4200)

# Get reads from a fastQ file, whether gzipped or plain text
def readFastQ(path, limit=None):
    of = gzip.open if path.endswith(".gz") else open
    with of(path, "rt") as handle:
        count = 0    
        for record in SeqIO.parse(handle, "fastq"):
            count += 1
            if limit and count > limit:
                break
            yield record

def readFastQFiles(fastQPaths, limitPerFile=None):
    for fastQPath in fastQPaths:
        for read in readFastQ(fastQPath, limitPerFile):
            yield (fastQPath, read)

def readFastQDirs(fastQDirs):
    for fastQDir in fastQDirs:
        for fastQPath in getFastQFiles(fastQDir):
            for read in readFastQ(fastQPath):
                yield read

def getFastQFiles(base, subdir=None):
    fastQDir = os.path.join(base, subdir) if subdir else base
    for file in sorted(filter(lambda f: f.endswith(".fastq.gz"), os.listdir(fastQDir))):
        yield os.path.join(fastQDir, file)

def getAllFastQDirs(fastQBaseDir):
    for dir in sorted(os.listdir(fastQBaseDir)):
        if dir.startswith("barcode"):
            yield os.path.join(fastQBaseDir, dir)

# Given a dictionary of sample name / fastq directory pairs, return a stream of
# (sample name, fastq directory) pairs
def getAllSampleDirs(fastQBaseDirs):
    for (samplePrefix, fastQBaseDir) in fastQBaseDirs.items():
        if samplePrefix:
            samplePrefix += "-"
        for fastQDir in getAllFastQDirs(fastQBaseDir):
            yield (samplePrefix + os.path.basename(fastQDir), fastQDir)

# Return all raw reads in a sub-directory on their own
def getReads(fastQDir):
    for (fastQPath,_) in getFastQAndHitsFiles(fastQDir):
        for read in readFastQ(fastQPath):
            yield read

def getAllReads(fastQBaseDir):
    for fastQDir in getAllFastQDirs(fastQBaseDir):
        for read in getReads(fastQDir):
            yield read

@dataclass
class Hit():
    primer: SeqRecord
    start: int     # read index, 0-based
    end: int       # read index past the end of the alignment
    rev: bool
    mr: float

# Set this to print alignments for debugging
primer_hits_to_print = 0

# Expand an ambiguous DNA sequence into all possible sequences
def expandAmbiguity(dna):
    seqs = [dna]
    for i, c in enumerate(dna):
        if c in IUPACData.ambiguous_dna_values:
            newSeqs = []
            for s in seqs:
                for n in IUPACData.ambiguous_dna_values[c]:
                    newSeqs.append(s[:i] + n + s[i+1:])
            seqs = newSeqs
    return seqs

def readPrimers(path, display=False, addRc=True, addRandom=True, expandAmbiguous=True):
    primers = []
    if display:
        print("Reading primers: " + path)
    for primer in SeqIO.parse(path, "fasta"):
        # The PairwiseAligner can't handle ambiguity codes, so expand out to multiple records
        seqs = expandAmbiguity(primer.seq) if expandAmbiguous else [primer.seq]
        if len(seqs) > 1:
            for i, s in enumerate(seqs):
                sfx = "." + str(i)
                primers.append(SeqRecord(s, id=primer.id + sfx, name=primer.name + sfx, description=primer.description + sfx))
        else:
            primers.append(primer)
        if display:
            print("  " + primer.description, end="")
            if len(seqs) > 1:
                print(" (%d variations)" % len(seqs), end="")
            print()
    if display:
        print("Read %d primers" % len(primers))

    if addRandom:
        # Add a random primer as a negative control
        minLen = min(len(p.seq) for p in primers)
        seq = Seq("".join([rnd.choice("ACGT") for i in range(minLen)]))
        primers.append(SeqRecord(seq, id="random", name="random", description="Random control"))

    # Precompute reverse complements for a ~5% speedup over using 'strand'
    if addRc:
        for primer in primers:
            primer.rcSeq = primer.seq.reverse_complement()

    # Store primer indicies for efficient serialization
    for i, primer in enumerate(primers):
        primer.index = i
        if ' ' in primer.description:
            primer.baseName = primer.description[:primer.description.rindex(' ')]

    return primers

# Get the amount of overlap between two hits. 1.0 means perfect overlap, 0.0 means not overlapping.
def computeOverlap(hit1, hit2):
    # The overlap fraction is the distance between the leftmost end the rightmost start
    # divided by the length of the smaller primer. This will be 1.0 for a 
    o = (min(hit1.end,hit2.end)-max(hit1.start,hit2.start)) / \
        min(len(hit1.primer),len(hit2.primer))
    return o if o > 0 else 0

def getPrimerAligner():
    aligner = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=0, gap_score=-1)
    return aligner

# Given a sequencing read, compute and return a list of primer match hits
# If allowOverlaps is true then overlapping results will be retained as long as their primer ID
# isn't the same up until the first "-" character.
def computePrimerHits(read, primers, allowOverlaps=False, matchThreshold=MATCH_THRESHOLD, aligner=None):
    aligner = getPrimerAligner() if aligner is None else aligner
    hits = []
    global primer_hits_to_print
    for primer in primers:
        seqs = [primer.seq]
        if hasattr(primer,'rcSeq'):
            seqs.append(primer.rcSeq)
        for query in seqs:
            # Find all alignments with a score above the threshold
            seqToMatch = read.seq
            while True:
                # Scoring seems to be ~10x faster than finding alignments
                score = aligner.score(seqToMatch, query)
                mr = round(score / len(query), ndigits=2)
                if mr < matchThreshold:
                    break
                alignment = aligner.align(seqToMatch, query)[0]
                assert alignment.score == score
                if primer_hits_to_print > 0:
                    print("Match: %.2f %s%s" % (mr, primer.description, "" if query==primer.seq else " (rev)"))
                    print(alignment)
                    primer_hits_to_print -= 1
                hit = Hit(
                    primer=primer, 
                    start=int(alignment.coordinates[0][0]), 
                    end=int(alignment.coordinates[0][-1]),
                    rev=query != primer.seq,
                    mr=mr)
                assert hit.start < len(read.seq) and hit.end <= len(read.seq), \
                    "Start %d, End %d, Len %d\n%s" % \
                    (alignment.start, alignment.end, len(read.seq), alignment)
                hits.append(hit)

                # Mask out the matched region
                if seqToMatch==read.seq:
                    seqToMatch=MutableSeq(read.seq)
                seqToMatch[hit.start:hit.end] = "X" * (hit.end-hit.start)
    # Remove redundant hits that are lower scoring
    hits.sort(key=lambda h: h.mr, reverse=True)
    trimmedHits = []
    for i in range(len(hits)):
        redundant = False
        for j in range(i):
            o = computeOverlap(hits[i],hits[j])
            if o >= OVERLAP_THRESHOLD:
                if not allowOverlaps or hits[i].primer.id.split("-")[0]==hits[j].primer.id.split("-")[0]:
                    redundant = True
                    break
        if not redundant:
            trimmedHits.append(hits[i])
    trimmedHits.sort(key=lambda h: h.start)
    return trimmedHits

# Given an alignment, extend it to completely cover the query
# Really what we want is a "semi-global" alignment, but that's not supported by Biopython.
# This gets us more-or-less the same effect, mainly for visualizing alignments. 
def extendAlignment(alignment):
    # Only works for forward alignments for now
    assert alignment.coordinates[1][0] < alignment.coordinates[1][-1]

    # Extend the ends of alignment to cover as much of the query as possible
    queryHead = min(alignment.coordinates[1][0], alignment.coordinates[0][0])
    queryTail = min(
        len(alignment.sequences[1]) - alignment.coordinates[1][-1],
        len(alignment.sequences[0]) - alignment.coordinates[0][-1])
    alignment.coordinates[0][0] -= queryHead
    alignment.coordinates[1][0] -= queryHead
    alignment.coordinates[0][-1] += queryTail
    alignment.coordinates[1][-1] += queryTail

    # Extend the alignment with gaps to ensure the entire query is covered,
    # for the case where the query runs off the start/end of the target.
    if(alignment.coordinates[1][0] > 0):
        alignment.coordinates = np.c_[[0,0], alignment.coordinates]
    if(alignment.coordinates[1][-1] < len(alignment.sequences[1])):
        alignment.coordinates = np.c_[alignment.coordinates, [len(alignment.sequences[0]), len(alignment.sequences[1])]]

def serializeHit(hit):
    return [hit.primer.index, hit.start, hit.end, hit.rev, hit.mr]

def deserializeHit(hitBuf, primers):
    return Hit(primer=primers[hitBuf[0]], start=hitBuf[1], end=hitBuf[2], rev=hitBuf[3], mr=hitBuf[4])

def getFastQPaths(fastQDir):
    for file in sorted(filter(lambda f: f.endswith(".fastq.gz"), os.listdir(fastQDir))):
        yield os.path.join(fastQDir, file)

def getHitsForFastQ(fastQPath):
    if fastQPath.endswith(".fastq.gz"):
        base = fastQPath.removesuffix(".fastq.gz")
    elif fastQPath.endswith(".fastq"):
        base = fastQPath.removesuffix(".fastq")
    return base+"-hits.json"

def getFastQAndHitsFiles(fastQDir):
    for fastQPath in getFastQPaths(fastQDir):
        hitsPath = getHitsForFastQ(fastQPath)
        yield (fastQPath, hitsPath)

# Find primer matches and save them to files if files don't already exist
def generateHitsFile(primers, fastQDir, overwrite=False):
    for fastQPath in getFastQPaths(fastQDir):
        generateSingleHitsFile(primers, fastQPath, overwrite)

def generateSingleHitsFile(primers, fastQPath, overwrite=False, limit=None, minQual=None, maxLength=None):
    hitsPath = getHitsForFastQ(fastQPath)
    if overwrite or not os.path.exists(hitsPath):
        print("Processing ", os.path.basename(fastQPath), end="")
        reads = 0
        start = time.process_time()
        serializedHitsPerRead = {}
        for read in readFastQ(fastQPath):
            reads += 1
            if reads % 100 == 0:
                print(".",end="")
            if maxLength and len(read.seq) > maxLength:
                continue
            if minQual and np.mean(read.letter_annotations["phred_quality"]) < minQual:
                continue
            hits = computePrimerHits(read, primers)
            serializedHitsPerRead[read.id]=[serializeHit(hit) for hit in hits]
            if limit and reads >= limit:
                break
            
        elapsed = time.process_time() - start
        print("  %.2fs" % elapsed)

        with open(hitsPath, "w") as f:
            json.dump(serializedHitsPerRead, f)
    else:
        print("Found ", os.path.basename(hitsPath))

# Stream all reads for a given subdirectory, along with the pre-computed primer matches
def getPrimerMatches(primers, fastQDir):
    if os.path.isdir(fastQDir):
        files = getFastQFiles(fastQDir)
    else:
        files = [fastQDir]

    for fastQPath in files:
        hitsPath = getHitsForFastQ(fastQPath)
        with open(hitsPath, "r") as hitsFile:
            serializedHitsPerRead = json.load(hitsFile)
            for read in readFastQ(fastQPath):
                if read.id in serializedHitsPerRead:
                    hits = []
                    for hitBuf in serializedHitsPerRead[read.id]:
                        hit = deserializeHit(hitBuf, primers)
                        assert hit.start < len(read.seq)
                        assert hit.end <= len(read.seq)
                        hits.append(hit)
                    yield (read, hits)

# Given a list of primers and a dictionary of sample prefixes to directories,
# stream tuples of (sample, read, hit) values for all reads in all samples.
# If given just a single string for fastQBaseDirs, then process just that directory
# with no sample prefix.      
def getAllPrimerMatches(primers, fastQBaseDirs):
    if isinstance(fastQBaseDirs, str):
        fastQBaseDirs = {"": fastQBaseDirs}
    
    for (sample, fastQDir) in getAllSampleDirs(fastQBaseDirs):
        for (read ,hit) in getPrimerMatches(primers, fastQDir):
            yield (sample, read, hit)

# Actually generate all the hits files
def generateAllHitsFiles(primers, fastQBaseDir):
    for fastQDir in getAllFastQDirs(fastQBaseDir):
        generateHitsFile(primers, fastQDir)

# Stream all plausible primer pairs
# This may include multiple overlapping pairs for a given read
def getPrimerPairs(primers, fastQPath, subdir=None):
    # If fastQPath is a list of paths then use that directly
    if isinstance(fastQPath, list):
        gen = ((os.path.basename(path), read, hits) for path in fastQPath for (read, hits) in getPrimerMatches(primers, path)) 
    elif subdir:
        gen = ((subdir, read, hits) for (read, hits) in getPrimerMatches(primers, os.path.join(fastQPath, subdir)))
    else:
        gen = getAllPrimerMatches(primers, fastQPath)
    for (subdir, read, hits) in gen:
        for hit1 in hits:
            if not hit1.rev:
                desc = hit1.primer.description[:hit1.primer.description.rindex(' ')]
                for hit2 in hits:
                    if hit2.rev and hit1.primer.baseName == hit2.primer.baseName \
                            and hit1.primer.description != hit2.primer.description:
                        span = hit2.start - hit1.end - 1 if hit2.start > hit1.start else hit1.start - hit2.end - 1
                        (p1,p2) = (hit1.primer.description, hit2.primer.description)
                        if len(p1) > len(p2) or (len(p1) == len(p2) and p1 > p2):
                            (p1,p2) = (p2,p1)
                        pairname = p1 + "/" + p2[p2.rindex(" ")+1:]
                        yield (subdir, read, hit1, hit2, span, pairname)

@dataclass
class SeqHit():
    target: SeqRecord
    strand: str
    score: int
    readStart: int = -1     # read index, 0-based
    readEnd: int = -1       # read index past the end of the alignment
    targetStart: int = -1
    targetEnd: int = -1

SEQ_SCORE_THRESHOLD = 100

# Given a sequencing read, look for matches to any of the given sequences.
# Rerturns a list of SeqHits
def seqMatch(read, seqs):
    # Using a gap score of -2 is almost twice as slow as -1, but avoids a lot of false positives
    al = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=-1, gap_score=-2)
    hits = []
    for target in seqs:
        for strand in ("+", "-"):
            # Find the single best alignment
            # Scoring is about 10x faster than aligning
            score = al.score(read.seq, target.seq, strand)
            if score >= SEQ_SCORE_THRESHOLD:
                # Resolve the full alignment
                alignment = al.align(read.seq, target.seq, strand)[0]
                assert alignment.score == score
                hit = SeqHit(
                    target=target, 
                    score=score,
                    strand=strand,
                    readStart = int(alignment.coordinates[0][0]),
                    readEnd = int(alignment.coordinates[0][-1]),
                    targetStart = int(alignment.coordinates[1][0]),
                    targetEnd = int(alignment.coordinates[1][-1])
                    )
                assert hit.readStart < len(read.seq) and hit.readEnd <= len(read.seq), \
                    "Start %d, End %d, Len %d\n%s" % \
                    (alignment.readStart, alignment.readEnd, len(read.seq), alignment)
                hits.append(hit)

    hits.sort(key=lambda h: h.score, reverse=True)
    return hits
