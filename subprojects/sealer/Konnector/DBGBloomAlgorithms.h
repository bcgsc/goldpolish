/**
 * Algorithms for a de Bruijn Graph using a Bloom filter
 * Copyright 2014 Canada's Michael Smith Genome Science Centre
 */
#ifndef DBGBLOOMALGORITHMS_H
#define DBGBLOOMALGORITHMS_H 1

#include "Common/Kmer.h"
#include "Common/KmerIterator.h"
#include "DBGBloom.h"
#include "Common/StringUtil.h"
#include "Common/Sequence.h"
#include "DataLayer/FastaReader.h"
#include "Graph/Path.h"
#include "btllib/nthash.hpp"
#include <climits>
#include <string>
#include <algorithm> // for std::max
#define NO_MATCH UINT_MAX

static inline Sequence pathToSeq(Path<Kmer> path)
{
	Sequence seq;
	assert(path.size() > 0);
	seq.append(path[0].str());
	for (unsigned i = 1; i < path.size(); i++)
		seq.append(1, path[i].getLastBaseChar());
	return seq;
}

/**
 * Choose a suitable starting kmer for a path search and
 * return its position. More specifically, find the kmer
 * closest to the end of the given sequence that is followed by
 * at least (numMatchesThreshold - 1) consecutive kmers that
 * are also present in the Bloom filter de Bruijn graph. If there
 * is no sequence of matches of length numMatchesThreshold,
 * use the longest sequence of matching kmers instead.
 *
 * @param seq sequence in which to find start kmer
 * @param k kmer size
 * @param g de Bruijn graph
 * @param numMatchesThreshold if we encounter a sequence
 * of numMatchesThreshold consecutive kmers in the Bloom filter,
 * choose the kmer at the beginning of that sequence
 * @param anchorToEnd if true, all k-mers from end of sequence
 * up to the chosen k-mer must be matches. (This option is used when
 * we wish to preserve the original sequences of the reads.)
 * @return position of chosen start kmer
 */
template<typename Graph>
static inline unsigned getStartKmerPos(const Sequence& seq,
	unsigned k, Direction dir, const Graph& g,
	unsigned numMatchesThreshold=1, bool anchorToEnd=false)
{
	assert(numMatchesThreshold > 0);

	if (seq.size() < k)
		return NO_MATCH;

	unsigned matchCount = 0;
	unsigned maxMatchLen = 0;
	unsigned maxMatchPos = 0;
	decltype(&btllib::NtHash::roll) roll_fn;
	size_t start_pos;
	int step;
	if (dir == FORWARD) {
		roll_fn = &btllib::NtHash::roll_back;
		start_pos = seq.length() - k;
		step = -1;
	} else {
		roll_fn = &btllib::NtHash::roll;
		start_pos = 0;
		step = 1;
	}
	int prev_pos = start_pos - step;
	btllib::NtHash nthash(seq, g.m_bloom.get_hash_num(), k, start_pos);
	while ((nthash.*roll_fn)()) {
		if ((int(nthash.get_pos()) - prev_pos > 1) || (!vertex_exists(nthash, g))) {
			if (matchCount > maxMatchLen) {
				assert(nthash.get_pos() - step >= 0 &&
					nthash.get_pos() - step < (int)(seq.length() - k + 1));
				maxMatchPos = nthash.get_pos() - step;
				maxMatchLen = matchCount;
			}
			if (anchorToEnd)
				break;
			matchCount = 0;			
		} else {
			matchCount++;
			if (matchCount >= numMatchesThreshold)
				return nthash.get_pos();
		}
		prev_pos = nthash.get_pos();
	}

	/* handle case where first/last kmer in seq is a match */
	if (matchCount > maxMatchLen) {
		assert(nthash.get_pos() >= 0 &&
			nthash.get_pos() < (int)(seq.length() - k + 1));
		maxMatchPos = nthash.get_pos();
		maxMatchLen = matchCount;
	}
	if (maxMatchLen == 0)
		return NO_MATCH;
	else
		return maxMatchPos;
}

struct BaseChangeScore {

	size_t m_pos;
	char m_base;
	unsigned m_score;

public:

	BaseChangeScore() :
			m_pos(0), m_base('N'), m_score(0){}

	BaseChangeScore(size_t pos, char base, unsigned score) :
			m_pos(pos), m_base(base), m_score(score){}

};

template<typename Graph>
static inline bool correctSingleBaseError(const Graph& g, unsigned k,
		FastaRecord& read, size_t& correctedPos, bool rc = false)
{
	if (read.seq.length() < k)
		return false;

	const std::string bases = "AGCT";
	const size_t minScore = 3;
	std::vector<BaseChangeScore> scores;

	for (size_t i = 0; i < read.seq.length(); i++) {

		size_t overlapStart = std::max((int) (i - k + 1), 0);
		size_t overlapEnd = std::min(i + k - 1, read.seq.length() - 1);
		assert(overlapStart < overlapEnd);
		Sequence overlapStr = read.seq.substr(overlapStart,
				overlapEnd - overlapStart + 1);
		size_t changePos = i - overlapStart;

		for (size_t j = 0; j < bases.size(); j++) {
			if (read.seq[i] == bases[j])
				continue;
			overlapStr[changePos] = bases[j];
			size_t score = 0;
			for (KmerIterator it(overlapStr, k, rc); it != KmerIterator::end();
					it++)
			{
				if (vertex_exists(*it, g))
					score++;
			}
			if (score > minScore)
				scores.push_back(BaseChangeScore(i, bases[j], score));
		}

	}

	if (scores.size() == 0)
		return false;

	BaseChangeScore bestScore;
	bool bestScoreSet = false;
	for (size_t i = 0; i < scores.size(); i++) {
		if (!bestScoreSet || scores[i].m_score > bestScore.m_score) {
			bestScore = scores[i];
			bestScoreSet = true;
		}
	}

	correctedPos = bestScore.m_pos;
	read.seq[correctedPos] = bestScore.m_base;

	return true;
}

/** Uppercase only bases that are present in original reads.
 *  @return number of mis-matching bases. */
static inline unsigned maskNew(const FastaRecord& read1,
		const FastaRecord& read2, FastaRecord& merged, int mask = 0)
{
	Sequence r1 = read1.seq, r2 = reverseComplement(read2.seq);
	if (mask) {
		transform(r1.begin(), r1.end(), r1.begin(), ::tolower);
		transform(r2.begin(), r2.end(), r2.begin(), ::tolower);
		transform(merged.seq.begin(), merged.seq.end(), merged.seq.begin(),
				::tolower);
	}
	unsigned mismatches = 0;
	for (unsigned i = 0; i < r1.size(); i++) {
		assert(i < merged.seq.size());
		if (r1[i] == merged.seq[i])
			merged.seq[i] = toupper(r1[i]);
		else
			mismatches++;
	}
	for (unsigned i = 0; i < r2.size(); i++) {
		assert(r2.size() <= merged.seq.size());
		unsigned merged_loc = i + merged.seq.size() - r2.size();
		if (r2[i] == merged.seq[merged_loc])
			merged.seq[merged_loc] = toupper(r2[i]);
		else
			mismatches++;
	}
	return mismatches;
}

#endif
