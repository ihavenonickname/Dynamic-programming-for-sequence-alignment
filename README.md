# Dynamic-programming-for-sequence-alignment

A toy code for sequences alignment, either local or global. If you have no idea, please read the following:
https://en.wikipedia.org/wiki/Sequence_alignment#Dynamic_programming

```
from seqalign import SequencesAligner

sa = SequencesAligner(seq1="GITHUB", seq2="GITTYHUB")
sa.initmatrix()
sa.mountmatrix()
print "\n".join(sa.alignment())
```
