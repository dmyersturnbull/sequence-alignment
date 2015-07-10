package com.github.dmyersturnbull.alignment;

import com.google.common.base.Preconditions;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

import javax.annotation.Nonnegative;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.Immutable;
import javax.validation.constraints.Size;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.ParameterizedType;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;

/**
 * Performs global and local alignment on DNA sequences.
 *
 * Example usage:
 * {@code
 *     SequenceAligner<DnaSequence, NucleotideCompound> aligner
 *     = new SequenceAligner<>(new SimpleGapPenalty(11, 1),
 *                                   SubstitutionMatrixHelper.getNuc4_4(),
 *                                   Alignments.PairwiseSequenceAlignerType.GLOBAL,
 *                                   DNASequence::new
 *                                  );
 *     SequenceAlignmentWithPvalue<DnaSequence, NucleotideCompound> alignment
 *     = aligner.alignAndCalcPvalue(200, sequenceA, sequenceB);
 * }
 *
 * @author Douglas Myers-Turnbull
 */
@Immutable
public class SequenceAligner<S extends AbstractSequence<C>, C extends Compound> {

	private static final int sf_defaultGapOpenPenalty = 11;
	private static final int sf_defaultGapExtensionPenalty = 1;
	private static final SubstitutionMatrix<NucleotideCompound> sf_defaultMatrix = SubstitutionMatrixHelper.getNuc4_4();

	private final GapPenalty m_gapPenalty;

	private final Function<String, S> m_creator;

	private final SubstitutionMatrix<C> m_matrix;
	private final Alignments.PairwiseSequenceAlignerType m_type;

	@Nonnull
	public static SequenceAligner<DNASequence, NucleotideCompound> withDefaultOptions(@Nonnull Alignments.PairwiseSequenceAlignerType type) {
		return new SequenceAligner<>(new SimpleGapPenalty(sf_defaultGapOpenPenalty, sf_defaultGapExtensionPenalty), sf_defaultMatrix, type, DNASequence::new);
	}

	public SequenceAligner(@Nonnull GapPenalty gapPenalty, @Nonnull SubstitutionMatrix<C> matrix,
	                       @Nonnull Alignments.PairwiseSequenceAlignerType type, @Nonnull SequenceCreator<S> creator) {
		m_gapPenalty = gapPenalty;
		m_matrix = matrix;
		m_type = type;
		m_creator = s -> {
			try {
				return creator.create(s);
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		};
	}


	public SequenceAlignmentWithPvalue<S, C> alignAndCalcPvalue(@Nonnegative int nSimulations, @Nonnull S a, @Nonnull S b) {
		SequenceAlignment<S, C> alignment = align(a, b);
		return calcPvalueByPermutation(nSimulations, alignment);
	}

	/**
	 * Permutes the target sequence.
	 * Note that this is only meaningful for sequences that are sufficiently diverse. For example, "TTTTTTT" against
	 * "TTTTTTT" will always receive a p-value of 0.
	 * To mitigate this issue, this method uses a slightly incorrect definition of p-value: the probability of
	 * getting an alignment score <em>more</em> extreme under the null. This isn't a problem for real applications.
	 */
	public SequenceAlignmentWithPvalue<S, C> calcPvalueByPermutation(@Nonnegative int nSimulations, @Nonnull SequenceAlignment<S, C> result) {
		//noinspection ConstantConditions
		Preconditions.checkNotNull(result.getSequencePair(), "SequenceAlignment result is null");
		int rank = 0;
		for (int i = 0; i < nSimulations; i++) {
			S permuted = m_creator.apply(
					randomlyPermute(result.getSequencePair().getOriginalSequences().get(1).toString())
			);
			// this is NOT strictly the definition of p-value, but it avoids issues with repetitive sequences
			if (result.getScore() > alignFast(result.getOriginalA(), permuted)) rank++;
		}
		return new SequenceAlignmentWithPvalue<>(result, 1d - 1d * rank / (nSimulations + 1d));
	}

	@Nonnull
	public SequenceAlignment<S, C> align(@Nonnull S a, @Nonnull S b) {
		PairwiseSequenceAligner<S, C> aligner = Alignments.getPairwiseAligner(a, b, m_type, m_gapPenalty, m_matrix);
		//noinspection ConstantConditions
		if (aligner == null) {
			throw new RuntimeException("Couldn't create aligner; something is wrong in Biojava");
		}
		// this actually performs the alignment; do it first
		SequencePair<S, C> pair = aligner.getPair();
		return new SequenceAlignment<>(toIntExact(aligner.getScore()), aligner.getSimilarity(), pair);
	}

	public int alignFast(@Nonnull S a, @Nonnull S b) {

		if (m_type == Alignments.PairwiseSequenceAlignerType.GLOBAL || m_type == Alignments.PairwiseSequenceAlignerType.GLOBAL_LINEAR_SPACE) {
			return calcAlignScoreFast(a, b, true);
		} else if (m_type == Alignments.PairwiseSequenceAlignerType.LOCAL || m_type == Alignments.PairwiseSequenceAlignerType.LOCAL_LINEAR_SPACE) {
			return calcAlignScoreFast(a, b, false);
		}
		throw new UnsupportedOperationException("Can't alignFast using type " + m_type);
	}

	private int calcAlignScoreFast(@Nonnull S a, @Nonnull S b, boolean global) {

		int bLength = b.getLength();
		int gop = m_gapPenalty.getOpenPenalty();
		int gep = m_gapPenalty.getExtensionPenalty();
		assert gop < 1 && gep < 1;

		int localBest = Integer.MIN_VALUE; // the best best (the best anywhere)
		int best = 0; // the best at the current cell; will be overwritten
		int[] currentM = new int[bLength + 1], currentX = new int[bLength + 1], currentY = new int[bLength + 1]; // will be overwritten
		int[] aboveM = new int[bLength + 1], aboveX = new int[bLength + 1], aboveY = new int[bLength + 1];
		aboveM[0] = m_matrix.getValue(a.getCompoundAt(1), b.getCompoundAt(1));
		aboveX[0] = Math.decrementExact(aboveM[0]);  // meaningless; just make sure it's not used and doesn't underflow
		for (int col = 1; col <= bLength; col++) {
			aboveX[col] = Math.addExact(m_gapPenalty.getOpenPenalty(), mulAndCheck(col, m_gapPenalty.getExtensionPenalty()));
			aboveM[col] = aboveY[col] = Math.decrementExact(min(aboveX[col], aboveX[col - 1])); // meaningless; just make sure they're not used and don't underflow
		}

		int leftY = gop + gep; // left of 1,1 (where we start) is 1,0, which must have been via Y
		int leftM = Math.decrementExact(leftY), leftX = Math.decrementExact(leftY);  // meaningless; just make sure it's not used and doesn't underflow

		for (int row = 1; row < a.getLength(); row++) {
			for (int col = 1; col < bLength; col++) {

				int match = m_matrix.getValue(a.getCompoundAt(row + 1), b.getCompoundAt(col + 1));

				currentY[col] = Math.addExact(
						gep,
						max(Math.addExact(gop, aboveM[col]), Math.addExact(gop, aboveX[col]), aboveY[col])
				);
				currentX[col] = Math.addExact(
						gep,
						max(Math.addExact(gop, leftM), Math.addExact(gop, leftY), leftX)
				);

				currentM[col] = Math.addExact(
						match,
						max(aboveM[col - 1], aboveX[col - 1], aboveY[col - 1])
				);

				if (!global && currentM[col] < 0) currentM[col] = 0;

				best = max(currentY[col], currentX[col], currentM[col]);
				if (best > localBest) localBest = best;

				leftM = currentM[col];
				leftX = currentX[col];
				leftY = currentY[col];
			}
			System.arraycopy(currentM, 0, aboveM, 0, aboveM.length);
			System.arraycopy(currentX, 0, aboveX, 0, aboveX.length);
			System.arraycopy(currentY, 0, aboveY, 0, aboveY.length);
		}

		if (global) return best;
		return localBest;
	}

	/**
	 * This only works if SequenceAligner is made abstract
	 */
	@SuppressWarnings("unchecked")
	private S create(String string) {
		try {
			return (S) ((Class<?>)((ParameterizedType)getClass().getGenericSuperclass())
					.getActualTypeArguments()[0])
					.getDeclaredConstructor(String.class)
					.newInstance(string);
		} catch (InstantiationException | IllegalAccessException | ClassCastException | NoSuchMethodException | InvocationTargetException e) {
			throw new RuntimeException("Couldn't create new instance", e);
		}
	}

	/**
	 * Converts {@code d} to an integer, and throws a {@link java.lang.ArithmeticException} if {@code d} is not an integer.
	 */
	private static int toIntExact(double d) {
		int x = (int)d;
		//noinspection FloatingPointEquality
		if (x != d) {
			throw new ArithmeticException("Cannot convert " + d + " to int");
		}
		return x;
	}

	@Nonnull
	private static String randomlyPermute(@Nonnull String string) {
		List<Character> list = new ArrayList<>();
		for (char c : string.toCharArray()) {
			list.add(c);
		}
		Collections.shuffle(list);
		StringBuilder sb = new StringBuilder();
		list.forEach(sb::append);
		return sb.toString();
	}

	/**
	 * @return The largest element in {@code array}
	 */
	private static int max(@Nonnull @Size(min = 1) int... array) {
		int max = Integer.MIN_VALUE;
		for (int x : array) {
			if (x > max) max = x;
		}
		return max;
	}

	/**
	 * @return The smallest element in {@code array}
	 */
	private static int min(@Nonnull @Size(min = 1) int... array) {
		int min = Integer.MAX_VALUE;
		for (int x : array) {
			if (x < min) min = x;
		}
		return min;
	}

	/**
	 * Multiply two integers, checking for overflow.
	 *
	 * Adapted from Apache Commons Math 3.6, {@code ArithmeticUtils.mulAndACheck(int, int)}.
	 *
	 * @param x Factor.
	 * @param y Factor.
	 * @return the product {@code x * y}.
	 * @throws ArithmeticException if the result can not be
	 * represented as an {@code int}.
	 */
	private static int mulAndCheck(int x, int y) throws ArithmeticException {
		long m = (long) x * (long) y;
		if (m >= -2147483648L && m <= 2147483647L) {
			return (int) m;
		} else {
			throw new ArithmeticException("Multiplying " + x + " with " + y + " would overflow");
		}
	}

	@FunctionalInterface
	public interface SequenceCreator<S> {
		S create(@Nonnull String string) throws Exception;
	}

}