#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <cassert>

#include "genome.h"
#include "genome1.h"

using namespace std;

int main(int argc, const char **argv)
{
 	if(argc == 1)
	{
		cout<<"usage: " << endl;
		cout<<"       " << argv[0] << " RPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " FPKM2TPM <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " shrink <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " filter <min-transcript-coverage> <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " top <integer-n> <in-gtf-file> <out-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " stats-exons <in-gtf-file> <exons-bins>"<<endl;
		cout<<"       " << argv[0] << " stats-length <in-gtf-file>"<<endl;
		cout<<"       " << argv[0] << " length <min-transcript-length> <max-transcript-length> <in-gtf-file> <out-gtf-file>"<<endl;
        cout<<"       " << argv[0] << " filter-tr-num <min-transcript-num> <in-gtf-file> <out-gtf-file>"<<endl;
        cout<<"       " << argv[0] << " select-tr <transcript-list> <in-gtf-file> <out-gtf-file>"<<endl;
        cout<<"       " << argv[0] << " TSSTES <in-gtf-file>"<<endl;
		cout<<"		  " << argv[0] << " remove-fp <in-gtf-file> <predictions-file> <out-gtf-file> <hard-mode>"<<endl;
		cout<<"		  " << argv[0] << " remove-fp-threshold <in-gtf-file> <predictions-file> <out-gtf-file> <hard-mode> <bp-threshold>"<<endl;
		cout<<"		  " << argv[0] << " filter-chrom <in-gtf-file> <chrom-file> <out-gtf-file>"<<endl;
		cout<<"		  " << argv[0] << " update-cov <in-gtf-file> <predictions-file> <out-gtf-file>"<<endl;
		cout<<"		  " << argv[0] << " update-tpm <in-gtf-file> <tpm-file> <out-gtf-file>"<<endl;
		cout<<"		  " << argv[0] << " get-coverage <in-gtf-file> <out-tsv-file>"<<endl;
		cout<<"		  " << argv[0] << " update-transcript-cov <in-gtf-file> <predictions-file> <out-gtf-file>"<<endl;
		return 0;
	}

	// cout << argv[1] << endl;
	
	if(string(argv[1]) == "FPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_FPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "RPKM2TPM")
	{
		genome gm(argv[2]);
		gm.assign_TPM_by_RPKM();
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "shrink")
	{
		genome1 gm;
		gm.shrink(argv[2]);
		gm.write(argv[3]);
	}

	if(string(argv[1]) == "filter")
	{
		genome1 gm;
		gm.filter(argv[3], atof(argv[2]));
		gm.write(argv[4]);
	}

	if(string(argv[1]) == "top")
	{
		genome1 gm;
		gm.top(argv[3], atoi(argv[2]));
		gm.write(argv[4]);
	}

	if(string(argv[1]) == "stats-exons")
	{
		genome1 gm;
		gm.stats_exons(argv[2], atoi(argv[3]));
	}

	if(string(argv[1]) == "stats-length")
	{
		genome1 gm;
		gm.stats_length(argv[2]);
	}

	if(string(argv[1]) == "length")
	{
		genome1 gm;
		gm.select_by_length(argv[4], atoi(argv[2]), atoi(argv[3]));
		gm.write(argv[5]);
	}
    
    if(string(argv[1]) == "filter-tr-num")
	{
		genome gm(argv[3]);
		set<string> geneIds = gm.filter_gene_with_minor_transcripts(atoi(argv[2]));
		gm.write_all(argv[3], argv[4], geneIds);
	}

    if(string(argv[1]) == "select-tr")
    {
        set<string> expressedTr;
        string tr;
        ifstream fin(argv[2]);//expressed list
        while(getline(fin, tr)) expressedTr.insert(tr);

        genome1 gm;
        gm.write_all(argv[3], argv[4], expressedTr);
    }

	if(string(argv[1]) == "TSSTES")
	{
		genome1 gm;
		gm.build_tsstes(argv[2]);
	}

	if(string(argv[1]) == "remove-fp" )
	{
		genome1 gm;
		assert(argc == 6);
		gm.filter_transcripts_with_tsstes(argv[2], argv[3], argv[4], atoi(argv[5]));
	}

	if(string(argv[1]) == "filter-chrom")
	{
		assert(argc == 5);
		genome1 gm;
		gm.filter_transcripts_by_chromosomes(argv[2], argv[3], argv[4]);
	}

	if(string(argv[1]) == "remove-fp-threshold" )
	{
		cout << "filtering based on threshold: " << atoi(argv[6]) << endl;
		genome1 gm;
		assert(argc == 7);
		gm.filter_transcripts_with_tsstes_threshold(argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]));
	}

	if (string(argv[1]) == "update-cov")
	{
		assert(argc == 5);
		genome1 gm;
		gm.update_coverage_by_prediction(argv[2], argv[3], argv[4]);
	}

	if (string(argv[1]) == "update-tpm")
	{
		assert(argc == 5);
		genome1 gm;
		gm.update_tpm(argv[2], argv[3], argv[4]);

	}

	if (string(argv[1]) == "get-cov")
	{
		assert(argc == 4);
		genome1 gm;
		gm.write_tss_tes_coverage(argv[2], argv[3]);
	}

	if (string(argv[1]) == "update-transcript-cov")
	{
		assert(argc == 5);
		// cout << "Updating transcript coverage based on TSS and TES predictions" << endl;
		genome1 gm;
		gm.update_transcript_coverage(argv[2], argv[3], argv[4]);
	}

	if (string(argv[1]) == "map-refseq-to-ucsc")
	{
		assert(argc == 5);
		// cout << "Mapping RefSeq to UCSC" << endl;
		genome1 gm;
		gm.map_refseq_to_ucsc(argv[2], argv[3], argv[4]);
	}

    return 0;
}
