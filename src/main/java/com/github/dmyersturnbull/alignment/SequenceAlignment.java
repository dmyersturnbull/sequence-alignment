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

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.Compound;

import javax.annotation.Nonnegative;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.Immutable;
import javax.validation.constraints.Max;
import javax.validation.constraints.Min;
import java.util.Objects;

/**
 * A result from sequence alignment.
 * @author Douglas Myers-Turnbull
 */
@Immutable
public class SequenceAlignment<S extends AbstractSequence<C>, C extends Compound> {

	private final int m_score;

	private final double m_similarity;
	private final SequencePair<S, C> m_sequencePair;

	public SequenceAlignment(int score, double similarity, @Nonnull SequencePair<S, C> sequencePair) {
		m_score = score;
		m_similarity = similarity;
		m_sequencePair = sequencePair;
	}

	public int getScore() {
		return m_score;
	}

	@Nonnegative
	public double getSimilarity() {
		return m_similarity;
	}

	@Nonnull
	public SequencePair<S, C> getSequencePair() {
		return m_sequencePair; // immutable; don't make a defensive copy
	}

	@Nonnegative
	public int getNMatches() {
		return m_sequencePair.getNumIdenticals();
	}

	@Nonnegative
	public int getNInsertionsInA() {
		int nOpens = 0;
		for (int i = 1; i <= m_sequencePair.getLength(); i++) {
			if (m_sequencePair.getCompoundInQueryAt(i).getShortName().equals("-")) nOpens++;
		}
		return nOpens;
	}

	@Min(0)
	@Max(1)
	public double getFractionIdenticalOfAligned() {
		int nId = 0;
		int nNotIndel = 0;
		for (int i = 1; i <= m_sequencePair.getLength(); i++) {
			String a = m_sequencePair.getCompoundInTargetAt(i).getShortName();
			String b = m_sequencePair.getCompoundInQueryAt(i).getShortName();
			if (!a.equals("-") && !b.equals("-")) {
				nNotIndel++;
				if (a.equals(b)) {
					nId++;
				}
			}
		}
		return (double)nId / nNotIndel;
	}

	@Nonnegative
	public int getNInsertionOpensInA() {
		int nOpens = 0;
		boolean insertionPrev = false;
		for (int i = 1; i <= m_sequencePair.getLength(); i++) {
			boolean isInsertion = m_sequencePair.getCompoundInQueryAt(i).getShortName().equals("-");
			if (isInsertion && !insertionPrev) nOpens++;
			insertionPrev = isInsertion;
		}
		return nOpens;
	}

	@Nonnegative
	public int getNInsertionOpensInB() {
		int nOpens = 0;
		boolean insertionPrev = false;
		for (int i = 1; i <= m_sequencePair.getLength(); i++) {
			boolean isInsertion = m_sequencePair.getCompoundInTargetAt(i).getShortName().equals("-");
			if (isInsertion && !insertionPrev) nOpens++;
			insertionPrev = isInsertion;
		}
		return nOpens;
	}

	@Nonnegative
	public int getNInsertionsInB() {
		int nOpens = 0;
		for (int i = 1; i <= m_sequencePair.getLength(); i++) {
			if (m_sequencePair.getCompoundInTargetAt(i).getShortName().equals("-")) nOpens++;
		}
		return nOpens;
	}

	@Nonnull
	public S getOriginalA() {
		return m_sequencePair.getOriginalSequences().get(0);
	}

	@Nonnull
	public S getOriginalB() {
		return m_sequencePair.getOriginalSequences().get(1);
	}

	@Nonnull
	public AlignedSequence<S, C> getAlignedA() {
		return m_sequencePair.getQuery();
	}

	@Nonnull
	public AlignedSequence<S, C> getAlignedB() {
		return m_sequencePair.getTarget();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("score=").append(m_score).append(System.lineSeparator());
		sb.append("nMatches=").append(getNMatches()).append(System.lineSeparator());
		sb.append("nInsertions=").append(getNInsertionsInA()).append(System.lineSeparator());
		sb.append("nDeletions=").append(getNInsertionsInB()).append(System.lineSeparator());
		sb.append("nInsertionOpens=").append(getNInsertionOpensInA()).append(System.lineSeparator());
		sb.append("nDeletionOpens=").append(getNInsertionOpensInB()).append(System.lineSeparator());
		sb.append("nDeletions=").append(getNInsertionsInB()).append(System.lineSeparator());
		sb.append("similarity=").append(String.format("%.6f", m_similarity)).append(System.lineSeparator());
		sb.append(getAlignedA()).append(System.lineSeparator());
		sb.append(getAlignedB()).append(System.lineSeparator());
		return sb.toString();
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) {
			return true;
		}
		if (o == null || getClass() != o.getClass()) {
			return false;
		}
		SequenceAlignment<?, ?> that = (SequenceAlignment<?, ?>) o;
		return Objects.equals(m_score, that.m_score)
				&& Objects.equals(m_similarity, that.m_similarity)
				&& Objects.equals(m_sequencePair, that.m_sequencePair);
	}

	@Override
	public int hashCode() {
		return Objects.hash(m_similarity, m_score, m_sequencePair.getQuery().toString(), m_sequencePair.getTarget().toString());
	}

}
