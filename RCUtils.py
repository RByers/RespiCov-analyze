from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Align
from Bio.Data import IUPACData
from dataclasses import dataclass

# Look for at least an 80% match against primers
# Below about 70% we seem to get huge numbers of matches just by chance
MATCH_THRESHOLD = 0.80

# Hits which overlap atleast 80% of a higher-score hit aren't reported
OVERLAP_THRESHOLD = 0.80

@dataclass
class Hit():
    primer: SeqRecord
    start: int     # read index, 0-based
    end: int       # read index past the end of the alignment
    rev: bool
    mr: float

# Set this to print alignments for debugging
primer_hits_to_print = 0

aligner = Align.PairwiseAligner(mode='local', match_score=1, mismatch_score=0, gap_score=-1)

def expandAmbiguity(dna):
    if len(dna) == 1:
        m = IUPACData.ambiguous_dna_values.get(dna)
        if m == None:
            return [dna]
        return [Seq(c) for c in m]
    head = expandAmbiguity(dna[:-1])
    tail = expandAmbiguity(dna[-1:])
    r = []
    for h in head:
        for t in tail:
            r.append(h+t)
    return r

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
        
    # Precompute reverse complements for a ~5% speedup over using 'strand'
    for primer in primers:
        primer.rcSeq = primer.seq.reverse_complement()
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
    hits = []
    global primer_hits_to_print
    for primer in primers:
        for query in (primer.seq, primer.rcSeq):
            # Find the single best alignment
            # Scoring seems to be ~10x faster than finding alignments
            score = aligner.score(read.seq, query)
            mr = round(score / len(query), ndigits=2)
            if mr >= MATCH_THRESHOLD:
                alignment = aligner.align(read.seq, query)[0]
                assert alignment.score == score
                # TODO consider looking for the same primer elsewhere by masking the match?
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
    return trimmedHits
