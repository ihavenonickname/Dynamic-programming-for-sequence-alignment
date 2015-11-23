# Dynamic-programming-for-sequence-alignment

Some basic stuff about dynamic programming. I've coded a python script with a class that gets two arrays and creates the alignment of them.

```
from seqalign import SequencesAligner

sa = SequencesAligner(seq1="GITHUB", seq2="GITTYHUB")
sa.initmatrix()
sa.mountmatrix()
print "\n".join(sa.alignment())
```
