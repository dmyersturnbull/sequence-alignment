# sequence-alignment

Simple Java code for sequence alignment and scoring, built on top of biojava-alignment.

In particular, can bootstrap to calculate a p-value by randomly permuting a sequence.
The bootstrapped alignments are calculated quickly using a sequence aligner that only records the alignment score (see `SequenceAligner.alignFast()`).
For example:

```java
SequenceAligner<DnaSequence, NucleotideCompound> aligner
    = new SequenceAligner<>(new SimpleGapPenalty(11, 1),
                            SubstitutionMatrixHelper.getNuc4_4(),
                            Alignments.PairwiseSequenceAlignerType.GLOBAL,
                            DNASequence::new
                           );
SequenceAlignmentWithPvalue<DnaSequence, NucleotideCompound> alignment
    = aligner.alignAndCalcPvalue(200, sequenceA, sequenceB);
```