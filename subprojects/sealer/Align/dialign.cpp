#include "dialign.h"
#include "Common/Options.h"
#include "Sequence.h"
#include <algorithm> // for min
#include <cassert>
#include <cmath> // for log
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime> // for clock
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

/** Score matrix. */
scr_matrix* smatrix;

/** Diagonal length probability distribution. */
prob_dist* pdist;

/** Return a DNA score matrix. */
static scr_matrix* newDefaultScoreMatrix()
{
	string s("ACGT?#$");
	struct scr_matrix* p = (scr_matrix*)calloc(1, sizeof *smatrix);
	p->length = s.size();
	p->num2char = (int*)calloc(256, sizeof(int));
	p->char2num = (int*)calloc(256, sizeof(int));
	for (unsigned i = 0; i < s.size(); ++i) {
		unsigned c = s[i];
		p->num2char[i] = c;
		p->char2num[c] = i;
	}

	p->data = (int*)calloc(s.size() * s.size(), sizeof(int));
	unsigned n = s.size() - 3; // ignore ?#$
	// Set the diagonal to 1.
	for (unsigned i = 0; i < n; i++)
		p->data[s.size() * i + i] = 1;
	p->max_score = 1;
	p->avg_sim_score = para->PROT_SIM_SCORE_THRESHOLD;
	p->dist = (int*)calloc(2, sizeof(int));
	p->dist[0] = n * n - n;
	p->dist[1] = n;
	return p;
}

/** Return a probability distribution for diagonal lengths
 * for a DNA score matrix.
 */
static prob_dist* newDefaultDiagProbDist()
{
	prob_dist *o = (prob_dist*)calloc(1, sizeof *o);
	o->smatrix = smatrix;
	unsigned length = 100;
	o->max_dlen = length;
	o->data = (long double**)calloc(
			length + 1, sizeof(long double *));
	o->log_data = (double**)calloc(length + 1, sizeof(double *));

	long double **dist = o->data;
	double **log_dist = o->log_data;
	const double* p = dna_diag_prob_100_exp_550000;
	for (unsigned i = 1; i <= length; i++) {
		unsigned mxscr = i * smatrix->max_score;
		dist[i] = (long double*)calloc(
				mxscr + 1, sizeof(long double));
		log_dist[i] = (double*)calloc(mxscr + 1, sizeof(double));
		for (unsigned scr = 0; scr <= mxscr; scr++) {
			double weight = *p++;
			assert(weight > 0);
			dist[i][scr] = weight;
			log_dist[i][scr] = -log(weight);
		}
	}
	return o;
}

/** Initialize dialign. */
void initDialign()
{
	// Score matrix
	smatrix = strlen(para->SCR_MATRIX_FILE_NAME) > 0
		? read_scr_matrix(para->SCR_MATRIX_FILE_NAME)
		: newDefaultScoreMatrix();
	if (para->DEBUG > 5)
		print_scr_matrix(smatrix);

	// Probability distribution for diagonal lengths
	pdist = strlen(para->DIAG_PROB_FILE_NAME) > 0
		? read_diag_prob_dist(smatrix, para->DIAG_PROB_FILE_NAME)
		: newDefaultDiagProbDist();
}

static void free_scr_matrix(struct scr_matrix* smatrix)
{
	free(smatrix->dist);
	free(smatrix->data);
	free(smatrix->char2num);
	free(smatrix->num2char);
	free(smatrix);
}

void free_prob_dist(struct prob_dist* pdist)
{
	unsigned int length = pdist->max_dlen;
	unsigned int i;
	for (i=1; i<=length; i++) {
		free(pdist->data[i]);
		free(pdist->log_data[i]);
	}
	free(pdist->data);
	free(pdist->log_data);
	free_scr_matrix(pdist->smatrix);
	free(pdist);
}

static void free_seq_col(struct seq_col* scol)
{
	unsigned int length = scol->length;
	unsigned int i;
	for (i=0; i<length; i++)
		free((scol->seqs[i]).data);
	free(scol->seqs);
	free(scol);
}

/** Print a dialign alignment. */
static ostream& print(ostream& out, const alignment& o,
		const string& consensus)
{
	const seq_col& scol = *o.scol;
	vector<int> proc(scol.length, 0);
	algn_pos **ap = o.algn;
	for (int s = 0; s < scol.length; s++) {
		const seq& sq = scol.seqs[s];
		for (int j = 0; j < o.max_pos; j++) {
			if (proc[s] < sq.length) {
				const algn_pos& ap1 = *find_eqc(ap, s, proc[s]);
				assert(j <= *ap1.eqcAlgnPos);
				if (*ap1.eqcAlgnPos == j) {
					char c = sq.data[proc[s]];
					if (toupper(c) == toupper(consensus[j]))
						out << '.';
					else if (ap1.state & para->STATE_ORPHANE)
						out << (char)tolower(c);
					else
						out << c;
					proc[s]++;
				} else
					out << '*';
			} else
				out << '*';
		}
		out << '\n';
	}
	return out;
}

/** Return the minimum number of matches. */
static unsigned countMatches(const alignment& o,
		const string& consensus)
{
	unsigned minMatches = consensus.size();
	const seq_col& scol = *o.scol;
	vector<int> proc(scol.length, 0);
	algn_pos **ap = o.algn;
	for (int s = 0; s < scol.length; s++) {
		unsigned matches = 0;
		const seq& sq = scol.seqs[s];
		for (int j = 0; j < o.max_pos; j++) {
			if (proc[s] < sq.length) {
				const algn_pos& ap1 = *find_eqc(ap, s, proc[s]);
				assert(j <= *ap1.eqcAlgnPos);
				if (*ap1.eqcAlgnPos == j) {
					char c = sq.data[proc[s]];
					if (toupper(c) == toupper(consensus[j]))
						matches++;
					proc[s]++;
				}
			}
		}
		minMatches = min(minMatches, matches);
	}
	return minMatches;
}

static struct seq_col* read_seqs(const vector<string>& amb_seqs)
{
	struct seq_col* scol = (struct seq_col*)calloc(1, sizeof(struct seq_col));
	struct seq* seqs = (scol->seqs = (struct seq*)calloc(amb_seqs.size(), sizeof(struct seq)));
	if(scol==NULL || seqs==NULL) {
		cerr << "read_seqs(): Out of memory !\n";
		exit(EXIT_FAILURE);
	}
	scol->length = amb_seqs.size();
	scol->avg_length = 0;

	seq* seq;
	for (size_t i=0; i<amb_seqs.size(); i++) {
		assert(!amb_seqs[i].empty());
		seq = &(scol->seqs[i]);
		seq->max_seen = 0;
		//seq->name = calloc(rlen, sizeof(char)); //do I need this?
		seq->num = i;
		seq->orf_frame=0;
		seq->crick_strand=0;
		//strncpy(seq->name, &(rline[1]), rlen-2);
		seq->data = (char*)calloc(amb_seqs[i].length()+1, sizeof(char));
		if (seq->data == NULL) {
			cerr << "seq->data out of memory !\n";
			exit(EXIT_FAILURE);
		}
		strcpy(seq->data, amb_seqs[i].c_str());
		seq->length = amb_seqs[i].length();
		scol->avg_length += amb_seqs[i].length();
		if(para->DEBUG >1) printf("DEBUG: seq:%s\n", seq->data);
	}
	scol->avg_length /= scol->length;
	if(para->DEBUG >1) printf("DEBUG: total # of amb_seqs: %i, avg_length: %i\n", scol->length, scol->avg_length);
	return scol;
}

// assume initial sequences contain only a/c/g/t/n
static string get_alignment_consensus(struct alignment *algn)
{
  struct seq_col *scol = algn->scol;
  unsigned int slen = scol->length;

  int j;
  unsigned int s,max;
  struct seq* sq;
  struct algn_pos **ap = algn->algn;

  prepare_alignment(algn);
  max = algn->max_pos;
  if (para->DEBUG > 5) printf("slen is %u, max pos is %u\n", slen, max);
  struct algn_pos *ap1;

	max = algn->max_pos;
	int* proc = new int[slen];
	for (j=0; j<(int)slen; j++)
		proc[j] = 0;
	string consensus;
	for (j=0; j<(int)max; j++) {
		char c = 'X';
		bool gap = false;
		for(s=0;s<slen;s++) {
			sq = &(scol->seqs[s]);
			if(proc[s] < sq->length) {
				ap1 = find_eqc(ap,s,proc[s]);
				if(*ap1->eqcAlgnPos==j) {
					char cur_char = toupper(sq->data[proc[s]]);
					c = c == 'X' ? cur_char
						: ambiguityOr(c, cur_char);
					proc[s]++;
				} else
					gap = true;
			} else
				gap = true;
		}
		consensus += gap ? tolower(c) : c;
	}
	delete[] proc;
	return consensus;
}

/** Align multiple sequences using DIALIGN-TX.
 * @param [out] alignment the alignment
 * @param [out] matches the minimum number of matches
 * @return the consensus sequence
 */
string dialign(const vector<string>& amb_seqs,
		string& alignment, unsigned& matches)
{
	int i;
	struct seq_col *in_seq_col = NULL;
	double tim = clock();

	in_seq_col = read_seqs(amb_seqs);

	// fast mode has higher threshold weights
	struct parameters *dialign_para = para;
	if(dialign_para->FAST_MODE)
		dialign_para->PROT_SIM_SCORE_THRESHOLD += 0.25;

	// Consider Anchors -> default for DNA: DO_ANCHOR = 0;
	struct alignment *algn = NULL;
	if (!dialign_para->FAST_MODE)
		algn = create_empty_alignment(in_seq_col);
	struct alignment *salgn = create_empty_alignment(in_seq_col);
	if (dialign_para->DEBUG > 1)
		printf("empty alignments created\n");

	// Compute pairwise diagonals
	struct diag_col *all_diags = find_all_diags(smatrix, pdist,
		in_seq_col, salgn, 1);
	double duration = (clock()-tim)/CLOCKS_PER_SEC;
	if (dialign_para->DEBUG > 1)
		printf("Found %i diags in %f secs\n",
			all_diags->diag_amount, duration);
	int diag_amount = all_diags->diag_amount;

	// Compute alignment
	double tim2 = clock();
	if (!dialign_para->FAST_MODE) {
		vector<diag*> cp_diags(all_diags->diag_amount);
		for(i = 0; i < diag_amount; i++) {
			cp_diags[i] = (diag*)malloc(sizeof(struct diag));
			*(cp_diags[i]) = *(all_diags->diags[i]);
		}
		guided_aligner(algn, in_seq_col, all_diags, smatrix,
			pdist, all_diags->gt_root, 1);

		for(i = 0; i < diag_amount; i++)
			all_diags->diags[i] = cp_diags[i];

		all_diags->diag_amount = diag_amount;
	}
	simple_aligner(in_seq_col, all_diags, smatrix, pdist,
		salgn, 1);
	duration = (clock()-tim2)/CLOCKS_PER_SEC;

	if (!dialign_para->FAST_MODE) {
		if (dialign_para->DEBUG > 1)
			printf("First alignment after %f secs. "
					"simple: %f guided: %f\n",
				duration, salgn->total_weight, algn->total_weight);
		else
			if (dialign_para->DEBUG > 1)
				printf("First alignment after %f secs. simple: %f \n",
					duration, salgn->total_weight);
	}

	free_diag_col(all_diags);

	dialign_para->DO_ANCHOR = 0; // anchors done

	// round 2+
	int round;
	char newFound = 0;
	int type;

	// consider sensitivity level
	if (!dialign_para->FAST_MODE) {
		if (dialign_para->SENS_MODE == 0) {
			dialign_para->DIAG_THRESHOLD_WEIGHT = 0.0;
		} else if (dialign_para->SENS_MODE == 1) {
			dialign_para->DIAG_THRESHOLD_WEIGHT
				= -log(0.75);//-log(.875+0.125/2.0);
		} else if (dialign_para->SENS_MODE == 2) {
			dialign_para->DIAG_THRESHOLD_WEIGHT
				= -log(0.5);//-log(0.875);
		}
	}

	int stype = (dialign_para->FAST_MODE ? 1 : 0);
	for (type = stype; type < 2; type++) {
		for (round = 2; round <= 20; round++) {
			tim2 = clock();
			all_diags = find_all_diags(smatrix, pdist,
				in_seq_col, (type ? salgn : algn), round);
			duration = (clock()-tim2)/CLOCKS_PER_SEC;
			if (dialign_para->DEBUG > 1)
				printf("Found %i diags after %f secs\n",
					all_diags->diag_amount, duration);
			if (all_diags->diag_amount == 0) {
				free_diag_col(all_diags);
				break;
			} else {
			// round 2 and further we use the simple aligner
				newFound = simple_aligner(in_seq_col,
					all_diags, smatrix, pdist,
					(type ? salgn : algn), round);
				free_diag_col(all_diags);
				if (!newFound)
					break;
			}
		}
	}
	if (dialign_para->DEBUG > 1)
		printf("Alignment ready!\n");

	if (!dialign_para->FAST_MODE) {
		if (dialign_para->DEBUG > 1)
			printf("Final alignment simple: %f guided: %f\n",
				salgn->total_weight, algn->total_weight);
	} else {
		if (dialign_para->DEBUG > 1)
			printf("Final alignment simple: %f \n",
				salgn->total_weight);
	}

	if (dialign_para->FAST_MODE
			|| salgn->total_weight > algn->total_weight) {
		if (!dialign_para->FAST_MODE)
			free_alignment(algn);
		algn = salgn;
	} else {
		free_alignment(salgn);
	}

	if (opt::verbose > 3)
		simple_print_alignment_default(algn);
	string consensus = get_alignment_consensus(algn);
	matches = countMatches(*algn, consensus);
	ostringstream ss;
	print(ss, *algn, consensus);
	alignment = ss.str();

	if (dialign_para->DEBUG > 0) {
		duration = (clock()-tim)/CLOCKS_PER_SEC;
		cerr << "Total time: " << duration << " s\n"
			"Total weight: " << algn->total_weight << '\n';
	}

	free_alignment(algn);
	free_seq_col(in_seq_col);
	return consensus;
}
