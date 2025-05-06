#include "genome1.h"
#include "config.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>

genome1::genome1()
{
}

int genome1::clear()
{
	transcripts.clear();
	intron_index.clear();
	return 0;
}

int genome1::build_all_transcripts(const string &file)
{
	genome gm(file);
	for (int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for (int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			transcripts.push_back(t);
		}
	}
	return 0;
}

int genome1::build_multiexon_transcripts(const string &file)
{
	genome gm(file);
	for (int i = 0; i < gm.genes.size(); i++)
	{
		const gene &g = gm.genes[i];
		for (int k = 0; k < g.transcripts.size(); k++)
		{
			const transcript &t = g.transcripts[k];
			if (t.exons.size() <= 1)
				continue;
			transcripts.push_back(t);
		}
	}
	return 0;
}

int genome1::build_intron_index()
{
	intron_index.clear();
	for (int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		assert(t.exons.size() >= 2);
		PI32 p = t.get_first_intron();
		if (intron_index.find(p) == intron_index.end())
		{
			set<int> s;
			s.insert(i);
			intron_index.insert(PPIS(p, s));
		}
		else
		{
			intron_index[p].insert(i);
		}
	}
	return 0;
}

int genome1::query(const transcript &t, const set<int> &fb)
{
	if (t.exons.size() <= 1)
		return -1;
	PI32 p = t.get_first_intron();
	if (intron_index.find(p) == intron_index.end())
		return -1;
	set<int> s = intron_index[p];
	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	for (int i = v.size() - 1; i >= 0; i--)
	{
		int k = v[i];
		transcript &x = transcripts[k];
		if (x.strand != t.strand)
			continue;
		if (x.seqname != t.seqname)
			continue;
		if (x.exons.size() != t.exons.size())
			continue;
		if (x.intron_chain_match(t) == false)
			continue;
		if (fb.find(k) != fb.end())
			continue;
		return k;
	}
	return -1;
}

int genome1::query(const transcript &t, int min_index)
{
	if (t.exons.size() <= 1)
		return -1;
	PI32 p = t.get_first_intron();
	if (intron_index.find(p) == intron_index.end())
		return -1;
	set<int> s = intron_index[p];
	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	for (int i = v.size() - 1; i >= 0; i--)
	{
		int k = v[i];
		if (k < min_index)
			continue;
		transcript &x = transcripts[k];
		if (x.strand != t.strand)
			continue;
		if (x.seqname != t.seqname)
			continue;
		if (x.exons.size() != t.exons.size())
			continue;
		if (x.intron_chain_match(t) == false)
			continue;
		return k;
	}
	return -1;
}

int genome1::shrink(const string &file)
{
	map<int, int> mm;
	vector<transcript> vv;
	build_multiexon_transcripts(file);
	int n = transcripts.size();
	if (n <= 0)
		return 0;

	sort(transcripts.begin(), transcripts.end(), transcript_cmp_intron_chain);

	vv.push_back(transcripts[0]);

	int c = 1;
	for (int k = 1; k < n; k++)
	{
		bool b = transcripts[k].intron_chain_match(transcripts[k - 1]);
		if (b == true)
			c++;
		if (b == true)
			continue;

		vv.push_back(transcripts[k]);
		if (mm.find(c) == mm.end())
			mm.insert(PII(c, 1));
		else
			mm[c]++;
		c = 1;
	}

	for (MII::iterator it = mm.begin(); it != mm.end(); it++)
	{
		printf("size %d : %d groups (by identifical intron-chain)\n", it->first, it->second);
	}

	transcripts = vv;
	return 0;
}

int genome1::select_by_length(const string &file, int min_length, int max_length)
{
	build_all_transcripts(file);
	vector<transcript> vv;
	for (int i = 0; i < transcripts.size(); i++)
	{
		if (transcripts[i].length() < min_length)
			continue;
		if (transcripts[i].length() > max_length)
			continue;
		vv.push_back(transcripts[i]);
	}
	transcripts = vv;
	return 0;
}

int genome1::filter(const string &file, double c)
{
	build_multiexon_transcripts(file);
	vector<transcript> vv;
	for (int i = 0; i < transcripts.size(); i++)
	{
		if (transcripts[i].coverage < c)
			continue;
		vv.push_back(transcripts[i]);
	}
	transcripts = vv;
	return 0;
}

int genome1::top(const string &file, int n)
{
	build_multiexon_transcripts(file);
	sort(transcripts.begin(), transcripts.end(), transcript_cmp_coverage);

	vector<transcript> vv;
	for (int i = transcripts.size() - 1; i >= 0; i--)
	{
		vv.push_back(transcripts[i]);
		if (vv.size() >= n)
			break;
	}
	transcripts = vv;
	return 0;
}

int genome1::compare(const genome1 &gy)
{
	MII x2y;
	MII y2x;
	set<int> fb;
	for (int i = gy.transcripts.size() - 1; i >= 0; i--)
	{
		const transcript &t = gy.transcripts[i];
		int k = query(t, fb);
		if (k == -1)
			continue;
		x2y.insert(PII(k, i));
		y2x.insert(PII(i, k));
		fb.insert(k);
	}

	int correct = x2y.size();
	int refsize = transcripts.size();
	int prdsize = gy.transcripts.size();

	double sen0 = correct * 100.0 / refsize;
	for (int i = 0; i < gy.transcripts.size(); i++)
	{
		double sen = correct * 100.0 / refsize;
		double pre = correct * 100.0 / (prdsize - i);

		if (sen * 2.0 < sen0)
			break;

		if (i % 100 == 0)
		{
			printf("ROC: reference = %d prediction = %d correct = %d sensitivity = %.2lf precision = %.2lf | coverage = %.3lf\n",
				   refsize, prdsize - i, correct, sen, pre, gy.transcripts[i].coverage);
		}

		if (y2x.find(i) != y2x.end())
			correct--;
	}
	return 0;
}

int genome1::stats_exons(const string &file, int n)
{
	build_all_transcripts(file);
	vector<int> counts;
	vector<int> length;
	counts.assign(n, 0);
	length.assign(n, 0);
	for (int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		int k = t.exons.size();
		int l = t.length();
		if (k >= n)
			k = n;
		counts[k - 1]++;
		length[k - 1] += l;
	}

	for (int k = 0; k < n; k++)
	{
		printf("transcripts with %d exons: %d\n", k + 1, counts[k]);
	}

	/*
	for(int k = 0; k < n; k++)
	{
		double ave = -1;
		if(counts[k] >= 1) ave = length[k] * 1.0 / counts[k];
	}

	for(int k = 0; k < n; k++)
	{
		if(counts[k] == 0) printf("0.0 ");
		else printf("%.2lf ", length[k] * 1.0 / counts[k]);
	}
	*/

	return 0;
}

int genome1::stats_length(const string &file)
{
	build_all_transcripts(file);
	map<int, int> m;
	for (int i = 0; i < transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		int l = t.length();

		printf("transcript has length of %d\n", l);
		if (m.find(l) == m.end())
			m.insert(pair<int, int>(l, 1));
		else
			m[l]++;
	}

	/*
	for(map<int,int>::iterator it = m.begin(); it != m.end(); it++)
	{
		printf("transcripts of lenght %d has %d\n", it->first, it->second);
	}
	*/

	return 0;
}

int genome1::print(int index)
{
	printf("genome %d: %lu transcripts, %lu distinct first intron\n", index, transcripts.size(), intron_index.size());
	return 0;
}

int genome1::write(const string &file)
{
	ofstream fout(file.c_str());
	for (int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].write(fout);
	}
	fout.close();
	return 0;
}

bool transcript_cmp_coverage(const transcript &x, const transcript &y)
{
	if (x.coverage < y.coverage)
		return true;
	else
		return false;
}

bool transcript_cmp_intron_chain(const transcript &x, const transcript &y)
{
	if (x.strand < y.strand)
		return true;
	if (x.strand > y.strand)
		return false;
	if (x.seqname < y.seqname)
		return true;
	if (x.seqname > y.seqname)
		return false;
	if (x.exons.size() < y.exons.size())
		return true;
	if (x.exons.size() > y.exons.size())
		return false;

	int n = x.exons.size() - 1;
	if (n >= 0 && x.exons[0].second < y.exons[0].second)
		return true;
	if (n >= 0 && x.exons[0].second > y.exons[0].second)
		return false;
	if (n >= 0 && x.exons[n].first < y.exons[n].first)
		return true;
	if (n >= 0 && x.exons[n].first > y.exons[n].first)
		return false;

	for (int k = 1; k < n - 1; k++)
	{
		int32_t x1 = x.exons[k].first;
		int32_t y1 = y.exons[k].first;
		int32_t x2 = x.exons[k].second;
		int32_t y2 = y.exons[k].second;
		if (x1 < y1)
			return true;
		if (x1 > y1)
			return false;
		if (x2 < y2)
			return true;
		if (x2 > y2)
			return false;
	}
	if (x.length() > y.length())
		return true;
	if (x.length() < y.length())
		return false;
	if (x.get_bounds().first < y.get_bounds().first)
		return true;
	if (x.get_bounds().first > y.get_bounds().first)
		return false;
	if (x.transcript_id.compare(y.transcript_id) < 0)
		return true;
	else
		return false;
}

string tostring(int32_t p)
{
	char buf[10240];
	sprintf(buf, "%d", p);
	return string(buf);
}

string compute_intron_hashing(const transcript &t)
{
	string h = "0";
	if (t.exons.size() <= 1)
		return h;
	int32_t p = t.exons[0].second;
	h.append(tostring(p));

	for (int k = 1; k < t.exons.size(); k++)
	{
		int32_t q1 = t.exons[k].first;
		int32_t q2 = t.exons[k].second;
		h.append(tostring(q1 - p));
		if (k == t.exons.size() - 1)
			break;
		h.append(tostring(q2 - q1));
		p = q2;
	}
	return h;
}

int genome1::write_all(const string &input, const string &output, set<string> expressedTr)
{
	if (input == "" || output == "")
		return 0;

	ifstream fin(input.c_str());
	if (fin.fail())
	{
		printf("open file %s error\n", input.c_str());
		return 0;
	}
	build_all_transcripts(input);

	ofstream fout(output.c_str());
	char line[102400];
	printf("Express %ld over %ld transcripts\n", expressedTr.size(), transcripts.size());
	while (fin.getline(line, 102400, '\n'))
	{
		item ge(line);
		if (expressedTr.find(ge.transcript_id) != expressedTr.end())
		{
			fout << line << endl;
		}
	}

	fin.close();
	fout.close();
	return 0;
}

int genome1::build_tsstes(const string &input)
{
	genome gm(input);

	map<pair<string, int32_t>, pair<int, int>> m1;
	map<pair<string, int32_t>, pair<int, int>> m2;
	for (int i = 0; i < gm.genes.size(); i++)
	{
		for (int j = 0; j < gm.genes[i].transcripts.size(); j++)
		{
			transcript &t = gm.genes[i].transcripts[j];
			// printf("TSS = %d, TES = %d, t-strand = %c\n", t.start, t.end, t.strand);

			pair<string, int32_t> p1 = make_pair(t.seqname, t.start + 1);
			pair<string, int32_t> p2 = make_pair(t.seqname, t.end);

			if (t.strand == '+')
			{
				if (m1.find(p1) == m1.end())
					m1.insert(make_pair(p1, pair<int, int>(1, 0)));
				else
					m1[p1].first++;
				if (m2.find(p2) == m2.end())
					m2.insert(make_pair(p2, pair<int, int>(1, 0)));
				else
					m2[p2].first++;
			}

			if (t.strand == '-')
			{
				if (m1.find(p2) == m1.end())
					m1.insert(make_pair(p2, pair<int, int>(0, 1)));
				else
					m1[p2].second++;
				if (m2.find(p1) == m2.end())
					m2.insert(make_pair(p1, pair<int, int>(0, 1)));
				else
					m2[p1].second++;
			}
		}
	}

	for (auto &x : m1)
		printf("TSS %s %d %d %d\n", x.first.first.c_str(), x.first.second, x.second.first, x.second.second);
	for (auto &x : m2)
		printf("TES %s %d %d %d\n", x.first.first.c_str(), x.first.second, x.second.first, x.second.second);

	return 0;
}

int genome1::filter_transcripts_with_tsstes(const string &input, const string &predictions, const string &output)
{
	build_all_transcripts(input);
	// read predictions tsv file
	std::unordered_map<std::string, std::unordered_set<int>> chrom_data_tss;
	std::unordered_map<std::string, std::unordered_set<int>> chrom_data_tes;
	std::vector<std::string> headers = {};

	std::ifstream pred_file(predictions);
	if (!pred_file)
	{
		std::cerr << "Error opening file: " << predictions << std::endl;
		return 1;
	}

	// Read headers
	std::string line, col_name;
	std::getline(pred_file, line); 
	std::stringstream ss_header(line);
	while(getline(ss_header,col_name, '\t')) headers.push_back(col_name);
	
	// Read data
	while (std::getline(pred_file, line))
	{
		std::stringstream ss(line);
		std::string value, chrom, site_type; 
		int pos, pred_value;
		int colIndex = 0;

		while (std::getline(ss, value, '\t') && colIndex < headers.size())
		{
			// columns[headers[colIndex]].push_back(value);
			switch (colIndex)
			{
			case 0:
				site_type = value;
				break;
			case 1: //chrom name
				chrom = value;
				break;
			case 2: // position
				pos = atoi(value.c_str());
				break;
			case 4: // prediction 
				pred_value = atoi(value.c_str());
				break;
			default:
				break;
			}
			colIndex++;
		}
		if (pred_value > 0)
		{
			if (site_type == "TSS")
			{
				if ( chrom_data_tss.find(chrom) == chrom_data_tss.end())
				{
					unordered_set<int> new_list;
					new_list.insert(pos);
					chrom_data_tss[chrom] = new_list;
				}
				else chrom_data_tss[chrom].insert(pos);
			}
			else if (site_type == "TES")
			{
				if ( chrom_data_tes.find(chrom) == chrom_data_tes.end())
				{
					unordered_set<int> new_list;
					new_list.insert(pos);
					chrom_data_tes[chrom] = new_list;
				}
				else chrom_data_tes[chrom].insert(pos);
			}
		}
	}
	for (const auto &x : chrom_data_tss)
	{
		cout << x.first << " ";
	}
	cout << endl;	
	for (const auto &y : chrom_data_tes)
	{
		cout << y.first << " ";
	}
	cout << endl;

	vector<transcript> filtered_tr;
	for( int i=0; i<(int)transcripts.size(); i++)
	{
		transcript &t = transcripts[i];
		if(chrom_data_tss.find(t.seqname) == chrom_data_tss.end() || chrom_data_tes.find(t.seqname) == chrom_data_tes.end()) 
		{
			cerr << "Can not find chromosome <" << t.seqname << "> in TSS or TES data" << endl;
			continue;
			//exit(1);
		}
		if(t.strand == '+')
		{
			if(chrom_data_tss[t.seqname].find(t.start+1) != chrom_data_tss[t.seqname].end() && chrom_data_tes[t.seqname].find(t.end) != chrom_data_tes[t.seqname].end())
			 filtered_tr.push_back(t);
		}
		else
		{
			if(chrom_data_tss[t.seqname].find(t.end) != chrom_data_tss[t.seqname].end() && chrom_data_tes[t.seqname].find(t.start+1) != chrom_data_tes[t.seqname].end())
			 filtered_tr.push_back(t);
		}
		
	} 
	transcripts = filtered_tr;
	write(output);
	return 0;
	
}


int genome1::filter_transcripts_by_chromosomes(const std::string &input_gtf, const std::string &chrom_file, const std::string &output_gtf)
{
    // Step 1: Load chromosome names from the file into a set
    std::unordered_set<std::string> allowed_chromosomes;
    std::ifstream chrom_ifs(chrom_file);
    if (!chrom_ifs)
    {
        std::cerr << "Error: cannot open chromosome list file: " << chrom_file << std::endl;
        return 1;
    }

    std::string chrom_line;
    while (std::getline(chrom_ifs, chrom_line))
    {
        if (!chrom_line.empty())
            allowed_chromosomes.insert(chrom_line);
    }
    chrom_ifs.close();

    // Step 2: Load all transcripts from the input GTF
    build_all_transcripts(input_gtf);

    // Step 3: Filter transcripts whose chromosome name is in the allowed list
    std::vector<transcript> filtered_transcripts;
    for (const auto &t : transcripts)
    {
        if (allowed_chromosomes.count(t.seqname))
            filtered_transcripts.push_back(t);
    }

    // Step 4: Write the filtered transcripts to the output GTF
    transcripts = filtered_transcripts;
    write(output_gtf);

    return 0;
}
