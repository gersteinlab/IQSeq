#include <assert.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>
#include "../common/read.h"

using namespace std;
using namespace boost;

int syntax_err() {
	L_(error) << "Usage:\ngroup_isos\n	log_level(0,1,2,...) gff_path gff_type out_format";
	return 1;
}

int lexical_cast_err() {
	L_(error) << "Lexical_cast error when converting arguments to numeric values";
	return 1;
}

int main(int argc, char * argv[])
{
	typedef tokenizer<char_separator<char> > tok;
	// initialize the logging system
	Log::ReportingLevel = LogLevel::info;
	
	// get parameters
	if (argc < 4) {
		return syntax_err();
	}
	int argi = 1;
	long log_level = lexical_cast<long>(argv[argi++]);
	Log::ReportingLevel = log_level;
	string gff_path = argv[argi++];
	string gff_type = argv[argi++];
	string out_format = argv[argi++];

	// global variables
	ifstream ifs;
	string line;

	set<string> inames;
	set<string> chrs;
	map<string, interval_list<long> > iname2il;
	map<string, string> iname2chr;
	map<string, string> iname2strand;

	map<string, set<string> > chrstrand2inames;
	map<string, interval_list<long> > chrstrand2groups;

	// load isoforms
	L_(info) << "Loading isoforms...";
	ifs.clear();
	ifs.open(gff_path.c_str());
	assert(ifs.is_open() && !ifs.eof());

	if (gff_type == "UCSC_GFF") {
		getline(ifs, line); // skip first line

		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			long start, end;
			string iname, tmp, chr, strand;
			iss >> chr >> tmp >> tmp >> start >> end
				>> tmp >> strand >> tmp >> iname;
			start--;
			inames.insert(iname);
			chrs.insert(chr);
			iname2chr[iname] = chr;
			iname2strand[iname] = strand;
			iname2il[iname].add_interval(start, end);
			chrstrand2inames[chr + strand].insert(iname);
		}
	} else if (gff_type == "GENELETS_GFF3") {
		getline(ifs, line); // skip first line
		getline(ifs, line); // skip second line
		while (getline(ifs, line) && !ifs.eof()) {
			istringstream iss(line);
			string type;
			long start, end;
			string infos;
			string chr, tmp, strand;
			iss >> chr >> tmp >> type;
			if (type == "exon") {
				iss >> start >> end
					>> tmp >> strand >> tmp >> infos;
				char_separator<char> sep(";");
				tok tInfo(infos, sep);
				BOOST_FOREACH (string const & info, tInfo) {
					if (info.size() > 7 &&
							info.substr(0, 7) == "Parent=") {
						string local_inames = info.substr(7);
						char_separator<char> sep2(",");
						tok tIname(local_inames, sep2);
						BOOST_FOREACH (string const & iname, tIname) {
							inames.insert(iname);
							chrs.insert("chr" + chr);
							iname2chr[iname] = "chr" + chr;
							iname2strand[iname] = strand;
							iname2il[iname].add_interval(start - 1, end);
							chrstrand2inames["chr"+chr+strand].insert(iname);
						}
					}
				}
			}
		}
	}
	L_(info) << "Loaded " << inames.size() << " isoforms...";

	BOOST_FOREACH(string const & iname, inames) {
		interval_list<long> const & il = iname2il[iname];
		long i_start = il.get_starts()[0];
		long i_end = il.get_ends()[il.get_num_intervals() - 1];
		chrstrand2groups[iname2chr[iname] + iname2strand[iname] ].
			add_interval(i_start, i_end);
	}

	vector<string> strands;
	strands.push_back("+");
	strands.push_back("-");
	unsigned long gene_count = 0;
	BOOST_FOREACH(string const & chr, chrs) {
		BOOST_FOREACH(string const & strand, strands) {
			interval_list<long> const & regions = chrstrand2groups[chr + strand];
			for (unsigned long i = 0;
					i < regions.get_num_intervals(); ++i) {
				long g_start = regions.get_starts()[i];
				long g_end = regions.get_ends()[i];
				BOOST_FOREACH(string const & iname, chrstrand2inames[chr + strand]) {
					interval_list<long> const & il = iname2il[iname];
					long start = il.get_starts()[0];
					long end = il.get_ends()[il.get_num_intervals() - 1];
					if (g_start <= start && end <= g_end) {
						if (out_format == "MAP") {
							cout << gene_count << "\t" << iname << endl;
						} else if (out_format == "GFF") {
							for (unsigned long j = 0; j < il.get_num_intervals(); ++j) {
								cout << iname2chr[iname] << "\t"
									<< "group_isos" << "\t"
									<< "exon" << "\t"
									<< (il.get_starts()[j] + 1) << "\t"
									<< il.get_ends()[j] << "\t"
									<< ".\t" << iname2strand[iname] << "\t"
									<< ".\t" << iname << endl;
							}
						}
					}
				}
				gene_count ++;
			}
		}
	}

	return 0;
}
