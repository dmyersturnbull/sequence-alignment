package com.github.dmyersturnbull.alignment;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SimpleSubstitutionMatrix;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.BeforeClass;
import org.junit.Test;

import java.net.URISyntaxException;

import static org.junit.Assert.assertEquals;

/**
 * A test for {@link SequenceAligner}.
 * @author Douglas Myers-Turnbull
 */
public class SequenceAlignerTest {

	private static SubstitutionMatrix<NucleotideCompound> sf_matrix;
	private static GapPenalty sf_gapPenalty;
	private static final int sf_m = 2;
	private static int sf_gop;
	private static int sf_gep;

	@BeforeClass
	public static void setUpBefore() throws URISyntaxException {
		sf_matrix = new SimpleSubstitutionMatrix<>(AmbiguityDNACompoundSet.getDNACompoundSet(), (short)sf_m, (short)-sf_m);
		sf_gapPenalty = new SimpleGapPenalty(5, 3); // and sf_m==2, so they're all coprime
		sf_gop = sf_gapPenalty.getOpenPenalty();
		sf_gep = sf_gapPenalty.getExtensionPenalty();
	}

	private static SequenceAligner<DNASequence, NucleotideCompound> getGlobalAligner() {
		return new SequenceAligner<>(sf_gapPenalty, sf_matrix, Alignments.PairwiseSequenceAlignerType.GLOBAL, DNASequence::new);
	}

	private static SequenceAligner<DNASequence, NucleotideCompound> getLocalAligner() {
		return new SequenceAligner<>(sf_gapPenalty, sf_matrix, Alignments.PairwiseSequenceAlignerType.LOCAL, DNASequence::new);
	}

	@Test
	public void testPerfect() throws Exception {
		DNASequence a = new DNASequence("ACTAACCGAGATTTTACCCCACGGTATTTTTT");
		DNASequence b = new DNASequence("ACTAACCGAGATTTTACCCCACGGTATTTTTT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		SequenceAlignment<DNASequence, NucleotideCompound> alignment = aligner.align(a, b);
		assertEquals(sf_m * "ACTAACCGAGATTTTACCCCACGGTATTTTTT".length(), alignment.getScore(), 0.00000001);
	}

	@Test
	public void testAlign() throws Exception {
		DNASequence a = new DNASequence("C GTAT  ATATCGCGCGC G CGATATATATATCT TCTCTAAAAAAA".replaceAll(" ", ""));
		DNASequence b = new DNASequence("G GTATATATATCGCGCGC A CGAT TATATATCTCTCTCTAAAAAAA".replaceAll(" ", ""));
//                                       --CGTATATATCGCGCGCGCGATATATATATCT-TCTCTAAAAAAA
//                                       GGTATATATATCGCGCGCACGAT-TATATATCTCTCTCTAAAAAAA
//  mismatches:                             ^              ^
// OR:
//		                                 CG--TATATATCGCGCGCGCGATATATATATCT-TCTCTAAAAAAA
//		                                 GGTATATATATCGCGCGCACGAT-TATATATCTCTCTCTAAAAAAA
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		SequenceAlignment<DNASequence, NucleotideCompound> alignment = aligner.align(a, b);
		int nMatches = "--CGTATATATCGCGCGCGCGATATATATATCT-TCTCTAAAAAAA".length() - 2 - 4;
		double expectedScore = nMatches * sf_m
				+				                       - 2 * sf_m // there are two mismatches
				+ + 3 * sf_gop + 4 * sf_gep; // there are 3 gap opens and either 1 or 4 extensions, depending on the def
		assertEquals(expectedScore, alignment.getScore(), 0.00000001);
		assertEquals(3, alignment.getNInsertionsInA());
		assertEquals(1, alignment.getNInsertionsInB());
		assertEquals(2, alignment.getNInsertionOpensInA());
		assertEquals(1, alignment.getNInsertionOpensInB());
	}

	//	@Test // TODO Biojava is even more wrong
	public void testAdjacentGaps() throws Exception {
		DNASequence a = new DNASequence("AAATTT  CCCATTT".replaceAll(" ", ""));
		DNASequence b = new DNASequence("AAA  GGGCCCATTT".replaceAll(" ", ""));
		GapPenalty gapPenalty = new SimpleGapPenalty(0, 0);
		@SuppressWarnings("unchecked")
		SequenceAligner<DNASequence, NucleotideCompound> aligner = new SequenceAligner<>(gapPenalty, sf_matrix, Alignments.PairwiseSequenceAlignerType.GLOBAL, DNASequence::new);
		int scoreA = aligner.align(a, b).getScore();
		assertEquals(sf_m * "AAA".length() + sf_m * "CCCATTT".length() + 2 * sf_gop - 5 * sf_gep, scoreA);
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "AAA".length() + sf_m * "CCCATTT".length() + 2 * sf_gop - 5 * sf_gep, score);
	}

	@Test
	public void testAlignFastPerfect() throws Exception {
		DNASequence a = new DNASequence("ACCNGGT");
		DNASequence b = new DNASequence("ACCNGGT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACCNGGT".length(), score);
	}

	@Test
	public void testAlignFastMismatches() throws Exception {
		DNASequence a = new DNASequence("ACTACT");
		DNASequence b = new DNASequence("ACGACG");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "AC".length() + sf_m * "AC".length() - 2 * sf_m, score);
	}

	@Test
	public void testAlignFastInsertion() throws Exception {
		DNASequence a = new DNASequence("ACTACTACTACTACT");
		DNASequence b = new DNASequence("ACTACTGACTACTACT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length() + sf_gop + sf_gep, score);
	}

	@Test
	public void testAlignFastDeletion() throws Exception {
		DNASequence a = new DNASequence("ACTACTGACTACTACT");
		DNASequence b = new DNASequence("ACTACTACTACTACT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length() + sf_gop + sf_gep, score);
	}

	@Test
	public void testAlignFastInsertionAndDeletion() throws Exception {
		DNASequence a = new DNASequence("ACTACTGACTACTACTGGTGGTGGGTGGAAAT   ".replaceAll(" ", ""));
		DNASequence b = new DNASequence("ACTACT ACTACTACTGGTGGTGG TGGAAATGGT".replaceAll(" ", ""));
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACT".length() + sf_m * "ACTACTACTGGTGGTGG".length() + sf_m * "TGGAAAT".length() +
				             3 * sf_gop + 5 * sf_gep, score);
	}

	@Test
	public void testAlignFastInsertionAtEnd() throws Exception {
		DNASequence a = new DNASequence("ACTACTACTACTACT");
		DNASequence b = new DNASequence("ACTACTACTACTACTGGGGGGGGG");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length() + sf_gop + "GGGGGGGGG".length() * sf_gep, score);
	}

	@Test
	public void testAlignFastDeletionAtEnd() throws Exception {
		DNASequence a = new DNASequence("ACTACTACTACTACTGGGGGGGGG");
		DNASequence b = new DNASequence("ACTACTACTACTACT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length() + sf_gop + "GGGGGGGGG".length() * sf_gep, score);
	}

	//	@Test // TODO
	public void testAlignFastDeletionAtStartLocal() throws Exception {
		DNASequence a = new DNASequence("GGGGGGGGGACTACTACTACTACT".trim());
		DNASequence b = new DNASequence("         ACTACTACTACTACT".trim());
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getLocalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length(), score);
	}

	@Test
	public void testAlignFastDeletionAtEndLocal() throws Exception {
		DNASequence a = new DNASequence("ACTACTACTACTACTGGGGGGGGG");
		DNASequence b = new DNASequence("ACTACTACTACTACT");
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getLocalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACTACTACT".length(), score);
	}

	@Test
	public void testAlignFastLocal() throws Exception {
		DNASequence a = new DNASequence("    ACTACTACTACTACTACT          ".replaceAll(" ", ""));
		DNASequence b = new DNASequence("GGGGACTACTGGGGGGGGGGGGGGGGGGGGGG".replaceAll(" ", ""));
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getLocalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACT".length(), score);
	}

	@Test
	public void testAlignFastLocalWithInsertion() throws Exception {
		DNASequence a = new DNASequence("    CCCCACTACTACT  ACTACTACTACTACTACTACT          ".replaceAll(" ", ""));
		DNASequence b = new DNASequence("GGGG    ACTACTACTGGACTACTACTGGGGGGGGGGGGGGGGGGGGGG".replaceAll(" ", ""));
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getLocalAligner();
		int score = aligner.alignFast(a, b);
		assertEquals(sf_m * "ACTACTACT".length() + sf_gop + 2 * sf_gep + sf_m * "ACTACTACT".length(), score);
	}

	//	@Test // TODO The score is one mismatch off
	public void testAlignFastHard() throws Exception {
		DNASequence a = new DNASequence("CGTAT  ATATCGCGCGCGCGATATATATATCT TCTCTAAAAAAA".replaceAll(" ", ""));
		DNASequence b = new DNASequence("GGTATATATATCGCGCGCACGAT TATATATCTCTCTCTAAAAAAA".replaceAll(" ", ""));
//		                                 CG--TATATATCGCGCGCGCGATATATATATCT-TCTCTAAAAAAA
//		                                 GGTATATATATCGCGCGCACGAT-TATATATCTCTCTCTAAAAAAA
		// mismatches:                   ^                 ^
		SequenceAligner<DNASequence, NucleotideCompound> aligner = getGlobalAligner();
		double expectedScore = sf_m * ("GTAT".length() + "ATATCGCGCGC".length() + "CGAT".length() +
				"TATATATCT".length() + "TCTCTAAAAAAA".length()) - 2 * sf_m + 3 * sf_gop + 4 * sf_gep;
		assertEquals(expectedScore, aligner.alignFast(a, b), 0);
	}

}