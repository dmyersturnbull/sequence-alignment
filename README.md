# sequence-alignment

A simple Java package for sequence alignment and scoring.

In particular, can calculate a p-value by a permutation test.
Alignments of permuted sequences are calculated quickly using an aligner that only records the score (see `SequenceAligner.alignFast()`).

Example usage:

```java
SequenceAligner<DnaSequence, NucleotideCompound> aligner =
	SequenceAligner.dna(Alignments.PairwiseSequenceAlignerType.GLOBAL)
	.setGapPenalty(new SimpleGapPenalty(11, 1))
	.setMatrix(SubstitutionMatrixHelper.getNuc4_4()) // the default anyway
	.build()
SequenceAlignmentWithPvalue<DnaSequence, NucleotideCompound> alignment
    = aligner.alignAndCalcPvalue(200, sequenceA, sequenceB);
    alignment.getPvalue(); // returns a double
```

**Warning: there is currently a bug in the p-value calculations; see [issue #1](https://github.com/dmyersturnbull/sequence-alignment/issues/1).**

**In addition, contiguous gaps are handled incorrectly; see [Biojava issue #213](https://github.com/biojava/biojava/issues/213).**

The software is licensed under the Apache License, Version 2.0.