#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include "../common/read.h"
#include "../common/fim.h"

using namespace std;
using namespace boost;

const string TASK_MLE = "mle";
const string TASK_FIM = "fim";
const string TYPE_KNOWN_ISO_ONLY = "KNOWN_ISO_ONLY";

bool matrices_are_equal(ublas::matrix<double> const & m1,
		ublas::matrix<double> const & m2) {
	unsigned long nrows = m1.size1();
	unsigned long ncols = m1.size2();
	if (nrows != m2.size1() || ncols != m2.size2()) {
		return false;
	}

	for (unsigned long p = 0; p < nrows - 1; ++p) {
		for (unsigned long q = 0; q < ncols - 1; ++q) {
			if (abs(m1(p,q) - m2(p,q)) > 1E-5) {
				return false;
			}
		}
	}
	return true;
}

int main(int argc, char * argv[])
{
	// initialize the logging system
	// Log::ReportingLevel = LogLevel::warning;
	Log::ReportingLevel = LogLevel::info;
	// Log::ReportingLevel = LogLevel::debug;
	// Log::ReportingLevel = LogLevel::debug2;
	
	if (argc < 2) {
		L_(info) << "usage:\nsim_iso task options";
		return 1;
	}
	
	int argi = 1;
	string task = argv[argi++];
	L_(info) << "Task: " << task;
	double unit_cost;
	unsigned long num_trials = 0;
	string isoform_path;
	double medium_read_cost;
	double short_read_cost;
	double short_pe_read_cost;
	double short_pe_insert_size_tolerance = 0;
	unsigned long min_partial_exon_size = 10;
	string read_type;
	bool known_iso_only = false;
	
	if (task == TASK_MLE) {
		// get parameters
		if (argc < 10) {
			L_(info) << "usage:\nsim_iso " << task << " unit_cost num_trials known_iso_only isoform_path medium_read_cost short_read_cost short_pe_read_cost short_pe_insert_size_tolerance [iso_probs]";
			return 1;
		}
		unit_cost = lexical_cast<double>(argv[argi++]);
		num_trials = lexical_cast<unsigned long>(argv[argi++]);
		known_iso_only = (argv[argi++] == TYPE_KNOWN_ISO_ONLY);
		isoform_path = argv[argi++];
		medium_read_cost = lexical_cast<double>(argv[argi++]);
		short_read_cost = lexical_cast<double>(argv[argi++]);
		short_pe_read_cost = lexical_cast<double>(argv[argi++]);
		short_pe_insert_size_tolerance = lexical_cast<double>(argv[argi++]);
	} else if (task == TASK_FIM) {
		if (argc < 6) {
			L_(info) << "usage:\nsim_iso " << task << " known_iso_only read_type isoform_path short_pe_insert_size_tolerance [iso_probs]";
			return 1;
		}
		known_iso_only = (argv[argi++] == TYPE_KNOWN_ISO_ONLY);
		read_type = argv[argi++];
		isoform_path = argv[argi++];
		short_pe_insert_size_tolerance = lexical_cast<double>(argv[argi++]);
	} else {
		L_(error) << "Unknown task: " << task;
		return 1;
	}

	// global variables
	ifstream ifs;
	const gsl_rng_type * rng_type;
	gsl_rng * rng;

	gsl_rng_env_setup();

	rng_type = gsl_rng_default;
	unsigned long seed = time(0) * getpid();
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set (rng, seed);

	// load insertion info
	ifs.clear();
	ifs.open(isoform_path.c_str());
	assert(ifs.is_open() && !ifs.eof());
	vec_gap iso_gaps = load_gaps_noheader(ifs, false);
	ifs.close();

	// generate exon set
	ExonSet exons;
	BOOST_FOREACH(ga_ptr const & gap, iso_gaps) {
		for (unsigned long i = 0; i < gap->exonCount; ++i) {
			exons.insert(gap->exonStarts[i], gap->exonEnds[i]);
		}
	}

	// build isoforms
	Isoforms iso;
	iso.build(exons, iso_gaps);
	iso.enumerate_possible_isoforms();
	iso.rebuild_possible_iso_probs();
	L_(info) << "Possible isoforms: " << iso.num_possible_isoforms
		<< ", exons: " << iso.num_exons
		<< ", known isoforms: " << iso.num_known_isoforms;
	L_(debug2) << "Known isoform indices: " << iso.known_iso_indices;
	L_(debug2) << "Exon lengths: " << iso.exon_lengths;
	L_(debug2) << "First known isoform exon indices: "
		<< iso.possible_iso_exon_indices[iso.known_iso_indices[0]];
	L_(debug2) << "First known isoform exon total lengths: "
		<< iso.possible_iso_exon_total_lengths[iso.known_iso_indices[0]];
	L_(debug2) << "First possible isoform exon indices: "
		<< iso.possible_iso_exon_indices[0];
	L_(debug2) << "First possible isoform exon total lengths: "
		<< iso.possible_iso_exon_total_lengths[0];
	L_(debug2) << "Second possible isoform exon indices: "
		<< iso.possible_iso_exon_indices[1];
	L_(debug2) << "Second possible isoform exon total lengths: "
		<< iso.possible_iso_exon_total_lengths[1];

	// set the iso.known_iso_probs
	unsigned int base_argi = argi;
	double total_prob = 0.0;
	while (argi < argc && argi - base_argi < iso.num_known_isoforms) {
		total_prob += lexical_cast<double>(argv[argi]);
		if (total_prob > 1.00001) {
			L_(warning) << "sum(input probs) > 1, only using the first "
				<< (argi - base_argi) << " values";
			break;
		}

		iso.known_iso_probs[argi - base_argi] =
			lexical_cast<double>(argv[argi]);
		argi++;
	}
	for (unsigned long i = argi - base_argi;
			i < iso.num_known_isoforms; ++i) {
		iso.known_iso_probs[i] = (1.0 - total_prob) /
			(double)(iso.num_known_isoforms - (argi - base_argi));
	}
	iso.rebuild_possible_iso_probs();

	/*
	// random
	shared_ptr<double> alpha(new double[iso.num_known_isoforms]);
	for (unsigned long i = 0; i < iso.num_known_isoforms; ++i) {
		(alpha.get())[i] = 0.5;
	}
	gsl_ran_dirichlet(rng, iso.num_known_isoforms,
			alpha.get(), iso.known_iso_probs); 
	iso.rebuild_possible_iso_probs();
	*/

	string probs = "";
	for (unsigned long i = 0; i < iso.num_known_isoforms; ++i) {
		probs += lexical_cast<string>(iso.known_iso_probs[i]) + ",";
	}
	L_(info) << "Known isoform probabilities: " << probs;

	// ARS
	shared_ptr<AccessibleReadStarts> mread_ars(
			new AccessibleReadStarts(
				iso.known_iso_exon_indices,
				iso.known_iso_exon_total_lengths,
				iso.exon_lengths,
				Read_stats::medium_single_expected_read_length));
	mread_ars->construct_ARS();
	shared_ptr<AccessibleReadStarts> sread_ars(
			new AccessibleShortReadStarts(
				iso.known_iso_exon_indices,
				iso.known_iso_exon_total_lengths,
				iso.exon_lengths,
				Read_stats::short_single_expected_read_length,
				min_partial_exon_size));
	sread_ars->construct_ARS();
	shared_ptr<AccessibleReadStarts> spread_ars(
			new AccessibleReadStarts(
				iso.known_iso_exon_indices,
				iso.known_iso_exon_total_lengths,
				iso.exon_lengths,
				Read_stats::short_single_expected_read_length * 2
				+ Read_stats::short_paired_expected_insert_size));
	spread_ars->construct_ARS();
	if (!known_iso_only) {
		mread_ars = shared_ptr<AccessibleReadStarts>(
				new AccessibleReadStarts(
					iso.possible_iso_exon_indices,
					iso.possible_iso_exon_total_lengths,
					iso.exon_lengths,
					Read_stats::medium_single_expected_read_length));
		mread_ars->construct_ARS();
		sread_ars = shared_ptr<AccessibleReadStarts>(
				new AccessibleShortReadStarts(
					iso.possible_iso_exon_indices,
					iso.possible_iso_exon_total_lengths,
					iso.exon_lengths,
					Read_stats::short_single_expected_read_length,
					min_partial_exon_size));
		sread_ars->construct_ARS();
		spread_ars = shared_ptr<AccessibleReadStarts>(
				new AccessibleReadStarts(
					iso.possible_iso_exon_indices,
					iso.possible_iso_exon_total_lengths,
					iso.exon_lengths,
					Read_stats::short_single_expected_read_length * 2
					+ Read_stats::short_paired_expected_insert_size));
		spread_ars->construct_ARS();
	}

	vector<double> iso_probs;
	unsigned long K;
	if (known_iso_only) {
		K = iso.num_known_isoforms;
		iso_probs.resize(K, 0);
		for (unsigned long k = 0; k < K; ++k) {
			iso_probs[k] = iso.known_iso_probs[k];
		}
	} else {
		K = iso.num_possible_isoforms;
		iso_probs = iso.possible_iso_probs;
	}

	if (task == TASK_MLE) {
		vector<vector<double> >  thetas;
		for (unsigned long t = 0; t < num_trials; ++t) {
			if ((t > 0) && (t % 10 == 0)) {
				L_(info) << "Executed " << t << " trials...";
			}
			vector<ublas::matrix<double> > m_delta_Gs;
			ublas::matrix<double> sum_fim = ublas::zero_matrix<double>(
					K - 1,
					K - 1);
			shared_ptr<Read> readp_medium_single(
					new Read_medium_single(mread_ars, rng));
			shared_ptr<Read> readp_short_single(
					new Read_short_single(sread_ars, rng));
			shared_ptr<Read> readp_short_paired(
					new Read_short_paired(spread_ars, rng,
						2 * Read_stats::short_single_expected_read_length,
						Read_stats::short_paired_expected_insert_size,
						short_pe_insert_size_tolerance));
			shared_ptr<Read> readp = readp_medium_single;
			for (int s = 0; s < 3; ++s) {
				// generate single reads
				unsigned long num_reads = 0;
				double ratio = 0;
				switch (s) {
					case 0:
						num_reads = (unsigned long)(
								medium_read_cost * unit_cost /
								Read_stats::medium_cost_per_bp /
								(double)Read_stats::medium_single_expected_read_length);
						readp = readp_medium_single;
						break;
					case 1:
						num_reads = (unsigned long)(
								short_read_cost * unit_cost /
								Read_stats::short_cost_per_bp /
								(double)Read_stats::short_single_expected_read_length);
						for (unsigned long i = 0;
								i < sread_ars->get_num_isoforms();
								++i) {
							if (known_iso_only) {
								ratio += iso.known_iso_probs[i]
									* (double)sread_ars->get_iso_ARS_total_length(i)
									/ (double)(iso.known_iso_exon_total_lengths[i][iso.known_iso_exon_indices[i].size() - 1] - Read_stats::short_single_expected_read_length);
							} else {
								ratio += iso.possible_iso_probs[i]
									* (double)sread_ars->get_iso_ARS_total_length(i)
									/ (double)(iso.possible_iso_exon_total_lengths[i][iso.possible_iso_exon_indices[i].size() - 1] - Read_stats::short_single_expected_read_length);
							}
						}
						num_reads = (unsigned long)(
								(double)num_reads * ratio);
						readp = readp_short_single;
						break;
					case 2:
						num_reads = (unsigned long)(
								short_pe_read_cost * unit_cost /
								Read_stats::short_cost_per_bp /
								(double)Read_stats::short_single_expected_read_length) / 2;
						readp = readp_short_paired;
						break;
					default:
						break;
				}
				ublas::matrix<double> m_delta_G = ublas::zero_matrix<double>(
						num_reads, K);
				generate_reads(num_reads, iso, readp, rng, m_delta_G, known_iso_only);
				L_(debug) << "Generated " << num_reads << " reads: " << readp->read_type;
				L_(debug2) << "Delta * G for the reads & possible isoforms : " << m_delta_G;
				m_delta_Gs.push_back(m_delta_G);
			}

			vector<double> theta;
			compute_iso_probs_EM(K, m_delta_Gs, theta);
			thetas.push_back(theta);
			L_(debug2) << "Estimated theta: " << theta;
		}
		L_(info) << "Executed " << num_trials << " trials... Done";

		double sum_var = 0;
		for (unsigned long k = 0; k < K - 1; ++k) {
			double sum = 0;
			double var = 0;
			for (unsigned long t = 0; t < num_trials; ++t) {
				sum += thetas[t][k];
				var += thetas[t][k] * thetas[t][k];
			}
			double mean = sum / num_trials;
			var /= (double)(num_trials - 1);
			var -= (double)(num_trials) / (double)(num_trials - 1) * mean * mean;
			sum_var += var;
		}

		double sum_var2 = 0;
		for (unsigned long k = 0; k < K - 1; ++k) {
			double var = 0;
			for (unsigned long t = 0; t < num_trials; ++t) {
				var += (thetas[t][k] - iso_probs[k]) *
					(thetas[t][k] - iso_probs[k]);
			}
			var /= (double)(num_trials);
			sum_var2 += var;
		}

		cout << medium_read_cost << "\t"
			<< short_read_cost << "\t"
			<< short_pe_read_cost << "\t"
			<< sqrt(sum_var / (double)(K - 1)) << "\t"
			<< sqrt(sum_var2 / (double)(K - 1))
			<< endl;
	} else if (task == TASK_FIM) {
		FIM fim;
		double ratio = 1;
		shared_ptr<Read> readp_medium_single(
				new Read_medium_single(mread_ars));
		shared_ptr<Read> readp_short_single(
				new Read_short_single(sread_ars));
		shared_ptr<Read> readp_short_pe(
				new Read_short_paired(spread_ars, rng,
					2 * Read_stats::short_single_expected_read_length,
					Read_stats::short_paired_expected_insert_size,
					short_pe_insert_size_tolerance));
		shared_ptr<Read> readp;
		if (read_type == "Read_medium_single") {
			readp = readp_medium_single;
		} else if (read_type == "Read_short_single") {
			readp = readp_short_single;
			ratio = 0;
			for (unsigned long i = 0;
					i < sread_ars->get_num_isoforms();
					++i) {
				if (known_iso_only) {
					ratio += iso.known_iso_probs[i]
						* (double)sread_ars->get_iso_ARS_total_length(i)
						/ (double)(iso.known_iso_exon_total_lengths[i][iso.known_iso_exon_indices[i].size() - 1] - Read_stats::short_single_expected_read_length);
				} else {
					ratio += iso.possible_iso_probs[i]
						* (double)sread_ars->get_iso_ARS_total_length(i)
						/ (double)(iso.possible_iso_exon_total_lengths[i][iso.possible_iso_exon_indices[i].size() - 1] - Read_stats::short_single_expected_read_length);
				}
			}
		} else if (read_type == "Read_short_paired") {
			readp = readp_short_pe;
		} else {
			L_(error) << "Unknown read type: " << read_type;
			return 1;
		}

		//		ublas::matrix<double> bf_fim_rst = fim.bruteforce_fim(
		//				readp, iso_probs);
		//		L_(info) << "BruteforceFIM(Read, PossibleIsoforms)";
		//		L_(debug) << bf_fim_rst;
		//		L_(info) << "# ofim calls: " << fim.get_ofim_call_count();

		for (int i = 0; i < 99; ++i) fim.fast_fim(readp, iso_probs);
		ublas::matrix<double> fast_fim_rst = fim.fast_fim(
				readp, iso_probs);
		L_(info) << "FastFIM(Read, PossibleIsoforms)";
		L_(debug)	<< fast_fim_rst;
		L_(info) << "# ofim calls: " << fim.get_ofim_call_count();
		//		if (matrices_are_equal(bf_fim_rst, fast_fim_rst)) {
		//			L_(info) << "BruteforceFIM = FastFIM";
		//		}

		//		for (int i = 0; i < 99; ++i) fim.faster_fim(readp, iso_probs);
		//		ublas::matrix<double> faster_fim_rst = fim.faster_fim(
		//				readp, iso_probs);
		//		L_(info) << "FasterFIM(Read, PossibleIsoforms)";
		//		L_(debug) << faster_fim_rst;
		//		L_(info) << "# ofim calls: " << fim.get_ofim_call_count();
		//		if (matrices_are_equal(fast_fim_rst, faster_fim_rst)) {
		//			L_(info) << "FastFIM = FasterFIM";
		//		}
		
		// output the fim for R to load
		// FIM::print_fim_for_R(cout, faster_fim_rst);
		FIM::print_fim_diag_for_R(cout, fast_fim_rst);
		FIM::print_fim_diag_for_R(cout, fast_fim_rst * ratio);
	} else {
		L_(error) << "Unknown task: " << task;
	}

	gsl_rng_free(rng);
	return 0;
}
