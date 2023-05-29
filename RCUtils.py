import os
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

# Get reads from a gzipped fastQ file
def readFastQ(path):
    with gzip.open(path, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            yield record

def getFastQFiles(base, subdir):
    fastQDir = os.path.join(base, subdir)
    for file in sorted(filter(lambda f: f.endswith(".fastq.gz"), os.listdir(fastQDir))):
        yield os.path.join(fastQDir, file)

def getAllFastQDirs(fastQBaseDir):
    for dir in sorted(os.listdir(fastQBaseDir)):
        if dir.startswith("barcode"):
            yield os.path.join(fastQBaseDir, dir)

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

def readPrimers(path):
    primers = []
    for primer in SeqIO.parse(path, "fasta"):
        # The PairwiseAligner can't handle ambiguity codes, so expand out to multiple records
        seqs = expandAmbiguity(primer.seq)
        if len(seqs) > 1:
            for i, s in enumerate(seqs):
                sfx = "." + str(i)
                primers.append(SeqRecord(s, id=primer.id + sfx, name=primer.name + sfx, description=primer.description + sfx))
        else:
            primers.append(primer)
        
    # Add a random primer as a negative control
    minLen = min(len(p.seq) for p in primers)
    seq = Seq("".join([random.choice("ACGT") for i in range(minLen)]))
    primers.append(SeqRecord(seq, id="random", name="random", description="Random control"))

    # Precompute reverse complements for a ~5% speedup over using 'strand'
    for primer in primers:
        primer.rcSeq = primer.seq.reverse_complement()

    # Store primer indicies for efficient serialization
    for i, primer in enumerate(primers):
        primer.index = i
        primer.baseName = primer.description[:primer.description.rindex(' ')]

    return primers

# Get the amount of overlap between two hits. 1.0 means perfect overlap, 0.0 means not overlapping.
def computeOverlap(hit1, hit2):
    # The overlap fraction is the distance between the leftmost end the rightmost start
    # divided by the length of the smaller primer. This will be 1.0 for a 
    o = (min(hit1.end,hit2.end)-max(hit1.start,hit2.start)) / \
        min(len(hit1.primer),len(hit2.primer))
    return o if o > 0 else 0

# Given a sequencing read, compute and return a list of primer match hits
# If allowOverlaps is true then overlapping results will be retained as long as their primer ID
# isn't the same up until the first "-" character.
def computePrimerHits(read, primers, allowOverlaps=False):
    aligner = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=0, gap_score=-1)
    hits = []
    global primer_hits_to_print
    for primer in primers:
        for query in (primer.seq, primer.rcSeq):
            # Find all alignments with a score above the threshold
            seqToMatch = read.seq
            while True:
                # Scoring seems to be ~10x faster than finding alignments
                score = aligner.score(seqToMatch, query)
                mr = round(score / len(query), ndigits=2)
                if mr < MATCH_THRESHOLD:
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
                seqToMatch[hit.start:hit.end] = "N" * (hit.end-hit.start)
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


def serializeHit(hit):
    return [hit.primer.index, hit.start, hit.end, hit.rev, hit.mr]

def deserializeHit(hitBuf, primers):
    return Hit(primer=primers[hitBuf[0]], start=hitBuf[1], end=hitBuf[2], rev=hitBuf[3], mr=hitBuf[4])

def getFastQAndHitsFiles(fastQDir):
    for file in sorted(filter(lambda f: f.endswith(".fastq.gz"), os.listdir(fastQDir))):
        fastQPath = os.path.join(fastQDir, file)
        yield (fastQPath, fastQPath.removesuffix(".fastq.gz")+"-hits.json")

# Find primer matches and save them to files if files don't already exist
def generateHitsFile(primers, fastQDir, overwrite=False):
    for (fastQPath, hitsPath) in getFastQAndHitsFiles(fastQDir):
        if overwrite or not os.path.exists(hitsPath):
            print("Processing ", os.path.basename(fastQPath), end="")
            reads = 0
            start = time.process_time()
            serializedHitsPerRead = []
            for read in readFastQ(fastQPath):
                reads += 1
                if reads % 100 == 0:
                    print(".",end="")
                hits = computePrimerHits(read, primers)
                serializedHitsPerRead.append([serializeHit(hit) for hit in hits])

            elapsed = time.process_time() - start
            print("  %.2fs" % elapsed)

            with open(hitsPath, "w") as f:
                json.dump(serializedHitsPerRead, f)
        else:
            print("Found ", os.path.basename(hitsPath))

# Stream all reads for a given subdirectory, along with the pre-computed primer matches
def getPrimerMatches(primers, fastQDir):
    for (fastQPath,hitsPath) in getFastQAndHitsFiles(fastQDir):
        with open(hitsPath, "r") as hitsFile:
            serializedHitsPerRead = json.load(hitsFile)
            
        for (readIdx,read) in enumerate(readFastQ(fastQPath)):
            hits = []
            for hitBuf in serializedHitsPerRead[readIdx]:
                hit = deserializeHit(hitBuf, primers)
                assert hit.start < len(read.seq)
                assert hit.end <= len(read.seq)
                hits.append(hit)
            yield (read, hits)
                
def getAllPrimerMatches(primers, fastQBaseDir):
    for fastQDir in getAllFastQDirs(fastQBaseDir):
        for (r,h) in getPrimerMatches(primers, fastQDir):
            yield (os.path.basename(fastQDir), r, h)

# Actually generate all the hits files
def generateAllHitsFiles(primers, fastQBaseDir):
    for fastQDir in getAllFastQDirs(fastQBaseDir):
        generateHitsFile(primers, fastQDir)

# Stream all plausible primer pairs
# This may include multiple overlapping pairs for a given read
def getPrimerPairs(primers, fastQBaseDir, subdir=None):
    if subdir:
        gen = ((subdir, read, hits) for (read, hits) in getPrimerMatches(primers, os.path.join(fastQBaseDir, subdir)))
    else:
        gen = getAllPrimerMatches(primers, fastQBaseDir)
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
class GenomeHit():
    target: SeqRecord
    strand: str
    score: int
    readStart: int = -1     # read index, 0-based
    readEnd: int = -1       # read index past the end of the alignment
    targetStart: int = -1
    targetEnd: int = -1

GENOME_SCORE_THRESHOLD = 100

# Given a sequencing read, look for matches to any of the given genomes.
# Rerturns either a single GenomeHit or None
def genomeMatch(read, genomes):
    # Using a gap score of -2 is almost twice as slow as -1, but avoids a lot of false positives
    al = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=-1, gap_score=-2)
    hits = []
    for target in genomes:
        for strand in ("+", "-"):
            # Find the single best alignment
            # Scoring is about 10x faster than aligning
            score = al.score(read.seq, target.seq, strand)
            if score >= GENOME_SCORE_THRESHOLD:
                # Resolve the full alignment
                alignment = al.align(read.seq, target.seq, strand)[0]
                assert alignment.score == score
                hit = GenomeHit(
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
