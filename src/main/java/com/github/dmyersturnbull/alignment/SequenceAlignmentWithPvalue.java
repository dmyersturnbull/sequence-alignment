/*
   Copyright 2015 Douglas Myers-Turnbull

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */

package com.github.dmyersturnbull.alignment;

import com.google.common.base.Objects;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

import javax.annotation.Nonnegative;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.Immutable;

/**
 * A sequence alignment result associated with a p-value.
 * @author Douglas Myers-Turnbull
 */
@Immutable
public class SequenceAlignmentWithPvalue<S extends AbstractSequence<C>, C extends Compound> extends SequenceAlignment<S, C> {

	private final double m_pvalue;

	public double getPvalue() {
		return m_pvalue;
	}

	public SequenceAlignmentWithPvalue(int score, @Nonnegative double similarity, @Nonnegative double pvalue, @Nonnull SequencePair<S, C> sequencePair) {
		super(score, similarity, sequencePair);
		m_pvalue = pvalue;
	}

	public SequenceAlignmentWithPvalue(@Nonnull SequenceAlignment<S, C> alignment, double pvalue) {
		this(alignment.getScore(), alignment.getSimilarity(), pvalue, alignment.getSequencePair());
	}

	@Override
	public String toString() {
		return super.toString() + System.lineSeparator() + "pvalue=" + m_pvalue;
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}
		//noinspection EqualsBetweenInconvertibleTypes
		if (!super.equals(o)) {
			return false;
		}
		SequenceAlignmentWithPvalue<?, ?> that = (SequenceAlignmentWithPvalue<?, ?>) o;
		return Objects.equal(m_pvalue, that.m_pvalue);
	}

	@Override
	public int hashCode() {
		return Objects.hashCode(super.hashCode(), m_pvalue);
	}
}
