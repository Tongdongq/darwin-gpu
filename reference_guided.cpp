

#include <iostream>
#include <string>
#include <cstring>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <map>
#include <thread>
#include <mutex>
#include "gact.h"
#include "fasta.h"
#include "ntcoding.h"
#include "seed_pos_table.h"
#include "ConfigFile.h"

#ifndef Z_COMPILE_USED
    #error "These files should be compiled using the z_compile.sh script"
#endif

#ifdef GPU
    #define BATCH 1
#endif

#define PRINT_SEED_POS_CONSTRUCT TRUE

//enum states {Z, D, I, M};

int NUM_BLOCKS;
int THREADS_PER_BLOCK;
int BATCH_SIZE;

// GACT scoring
int gact_sub_mat[10];
int gap_open;
int gap_extend;

// D-SOFT parameters
std::string seed_shape;
std::string seed_shape_str;
uint32_t bin_size;
int dsoft_threshold;
int num_seeds;
int seed_occurence_multiple;
int max_candidates;
int num_nz_bins;
bool ignore_lower;

// GACT first tile
int first_tile_size;
int first_tile_score_threshold;

//GACT extend 
int tile_size;
int tile_overlap;

//Multi-threading
int num_threads;

static std::string reference_string;
static std::string query_string;

uint32_t reference_length;
uint32_t query_length;

struct timeval start, end_time;

// for some reason, the reference_lenghts vector contains alternating values - zeroes
// same for reads_lenghts
std::vector<long long int> reference_lengths;
std::vector<std::string> reference_seqs;

std::vector<long long int> reads_lengths;
std::vector<std::string> reads_seqs;
std::vector<std::string> rev_reads_seqs;

std::vector<std::vector<std::string> > reference_descrips;
std::vector<long long int> reference_fileposs;

std::vector<std::vector<std::string> > reads_descrips;
std::vector<long long int> reads_fileposs;

char** reads_char;
char** rev_reads_char;

std::map<int, uint32_t> chr_id_to_start_bin;
std::map<uint32_t, int> bin_to_chr_id;

SeedPosTable *sa;

std::mutex io_lock;

std::map<char, char> rcmap;

void initRCMap() {
    rcmap['a'] = 'T';
    rcmap['A'] = 'T';
    rcmap['t'] = 'A';
    rcmap['T'] = 'A';
    rcmap['c'] = 'G';
    rcmap['C'] = 'G';
    rcmap['g'] = 'C';
    rcmap['G'] = 'C';
    rcmap['n'] = 'n';
    rcmap['N'] = 'N';
}

std::string RevComp(std::string seq) {
    std::string rc = "";
    for (int i = seq.size()-1; i >= 0; i--) {
        if (seq[i] != 'a' && seq[i] != 'A' &&
                seq[i] != 'c' && seq[i] != 'C' &&
                seq[i] != 'g' && seq[i] != 'G' &&
                seq[i] != 't' && seq[i] != 'T' &&
                seq[i] != 'n' && seq[i] != 'N') {
            std::cout<<"Bad Nt char: "<< seq[i] <<std::endl;
            exit(1);
        }
        rc += rcmap[seq[i]];
    }
    return rc;
}


void PrintTileLocation (std::string read_name, \
    uint32_t candidate_hit, uint32_t last_hit_offset, char strand) {

    int         chr_id    = bin_to_chr_id[candidate_hit/bin_size];
    std::string chr       = reference_descrips[chr_id][0];
    uint32_t    start_bin = chr_id_to_start_bin[chr_id];

    std::cout << \
    read_name << " " << \
    chr << " " << \
    candidate_hit - (start_bin*bin_size) << " " << \
    last_hit_offset << " " << \
    strand << std::endl;
}

// each CPU thread aligns a batch of multiple reads against the SeedTable
// the params start_read_num and last_read_num indicate which readIDs are to be aligned for each CPU thread
void AlignReads (int start_read_num, int last_read_num) {

    uint32_t log_bin_size = (uint32_t) (log2(bin_size));
    int num_bins = 1 + (reference_length >> log_bin_size);

    // candidate_hit_offset contains 64 bits
    // the high 32 bits contain the position in the ref, the low 32 bits the position in the read/query
    uint64_t* candidate_hit_offset;

    struct timeval begin, finish;

    //    gettimeofday(&begin, NULL);
    uint32_t* nz_bins_array = new uint32_t[num_nz_bins];
    uint64_t* bin_count_offset_array = new uint64_t[num_bins];
    candidate_hit_offset = new uint64_t[max_candidates];

    for(int i=0; i < num_bins; i++) {
        bin_count_offset_array[i] = 0;
    }

    //    gettimeofday(&finish, NULL);
    //    long useconds = finish.tv_usec - begin.tv_usec;
    //    long seconds = finish.tv_sec - begin.tv_sec;
    //    long mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;

#ifdef BATCH
    GACT_call *GACT_calls_for, *GACT_calls_rev;
#endif

    int num_candidates_for, num_candidates_rev;

    for (int k = start_read_num; k < last_read_num; k++) {
        int len = reads_lengths[k];

        // Forward reads
        num_candidates_for = sa->DSOFT(reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);
#ifdef BATCH
        GACT_calls_for = new GACT_call[num_candidates_for];
#endif
        io_lock.lock();
        std::cout << "Read (+) " << k << ": " << num_candidates_for << std::endl;
        for (int i = 0; i < num_candidates_for; i++) {
            PrintTileLocation(reads_descrips[k][0], \
                (candidate_hit_offset[i] >> 32), \
                ((candidate_hit_offset[i] << 32) >> 32), '+');
            int ref_pos = (candidate_hit_offset[i] >> 32);
            int chr_id = bin_to_chr_id[ref_pos/bin_size];
            uint32_t start_bin = chr_id_to_start_bin[chr_id];
            ref_pos -= start_bin*bin_size;
            int query_pos = candidate_hit_offset[i] & 0xffffffff;

#ifdef BATCH
            // store GACT_call
            GACT_calls_for[i].ref_id = chr_id;
            GACT_calls_for[i].query_id = k;
            GACT_calls_for[i].ref_pos = ref_pos;
            GACT_calls_for[i].query_pos = query_pos;
            GACT_calls_for[i].score = 0;
            GACT_calls_for[i].first = 1;
            GACT_calls_for[i].reverse = 1;
            GACT_calls_for[i].terminate = 0;
#else   // perform GACT immediately
            GACT((char*)reference_seqs[chr_id].c_str(), reads_char[k], \
                reference_lengths[chr_id], len, \
                tile_size, tile_overlap, \
                ref_pos, query_pos, first_tile_score_threshold);
#endif
        }   // end for all num_candidates_for seed hits
        io_lock.unlock();
#ifdef BATCH
        for(int i = 0; i < num_candidates_for; ++i){
        	GACT_call *c = &(GACT_calls_for[i]);
        	printf("GACT_call %d, ref_id: %d, query_id: %d, ref_pos: %d, query_pos: %d\n", i, c->ref_id, c->query_id, c->ref_pos, c->query_pos);
        }
#endif


        // Reverse complement reads
        int num_candidates_rev = sa->DSOFT(rev_reads_char[k], len, num_seeds, dsoft_threshold, candidate_hit_offset, bin_count_offset_array, nz_bins_array, max_candidates);
#ifdef BATCH
        GACT_calls_rev = new GACT_call[num_candidates_rev];
#endif
        io_lock.lock();
        std::cout << "Read (-) " << k << ": " << num_candidates_rev << std::endl;
        for (int i = 0; i < num_candidates_rev; i++) {
            PrintTileLocation(reads_descrips[k][0], \
                (candidate_hit_offset[i] >> 32), \
                ((candidate_hit_offset[i] << 32) >> 32), '-');
            int ref_pos = (candidate_hit_offset[i] >> 32);
            int chr_id = bin_to_chr_id[ref_pos/bin_size];
            uint32_t start_bin = chr_id_to_start_bin[chr_id];
            ref_pos -= start_bin*bin_size;
            int query_pos = candidate_hit_offset[i] & 0xffffffff;

#ifdef BATCH
            GACT_calls_rev[i].ref_id = chr_id;
            GACT_calls_rev[i].query_id = k;
            GACT_calls_rev[i].ref_pos = ref_pos;
            GACT_calls_rev[i].query_pos = query_pos;
            GACT_calls_rev[i].score = 0;
            GACT_calls_rev[i].first = 1;
            GACT_calls_rev[i].reverse = 1;
            GACT_calls_rev[i].terminate = 0;
#else   // perform GACT immediately
            GACT((char*)reference_seqs[chr_id].c_str(), reads_char[k], \
                reference_lengths[chr_id], len, \
                tile_size, tile_overlap, \
                ref_pos, query_pos, first_tile_score_threshold);
#endif
        }   // end for all num_candidates_rev seed hits
        io_lock.unlock();
    }   // end for every read assigned to this CPU thread

    printf("num_candidates: %d %d\n", num_candidates_for, num_candidates_rev);

#ifdef BATCH

    GACT_Batch(GACT_calls_for, num_candidates_for);

    delete[] GACT_calls_for;
    delete[] GACT_calls_rev;
#endif

    delete[] bin_count_offset_array;
    delete[] nz_bins_array;
}

int main(int argc, char *argv[]) {

    if (argc < 3) {
        fprintf(stderr, "Usage: ./reference_guided <REFERENCE>.fasta <READS>.fasta [NUM_BLOCKS THREADS_PER_BLOCK]\n");
        exit(1);
    }

    ConfigFile cfg("params.cfg");

    // GACT scoring
    gact_sub_mat[0] = cfg.Value("GACT_scoring", "sub_AA");
    gact_sub_mat[1] = cfg.Value("GACT_scoring", "sub_AC");
    gact_sub_mat[2] = cfg.Value("GACT_scoring", "sub_AG");
    gact_sub_mat[3] = cfg.Value("GACT_scoring", "sub_AT");
    gact_sub_mat[4] = cfg.Value("GACT_scoring", "sub_CC");
    gact_sub_mat[5] = cfg.Value("GACT_scoring", "sub_CG");
    gact_sub_mat[6] = cfg.Value("GACT_scoring", "sub_CT");
    gact_sub_mat[7] = cfg.Value("GACT_scoring", "sub_GG");
    gact_sub_mat[8] = cfg.Value("GACT_scoring", "sub_GT");
    gact_sub_mat[9] = cfg.Value("GACT_scoring", "sub_TT");
    gap_open        = cfg.Value("GACT_scoring", "gap_open");
    gap_extend      = cfg.Value("GACT_scoring", "gap_extend");

    // D-SOFT parameters
    seed_shape_str          = (std::string) cfg.Value("DSOFT_params", "seed_shape");
    bin_size                = cfg.Value("DSOFT_params", "bin_size");
    dsoft_threshold         = cfg.Value("DSOFT_params", "threshold");
    num_seeds               = cfg.Value("DSOFT_params", "num_seeds");
    seed_occurence_multiple = cfg.Value("DSOFT_params", "seed_occurence_multiple");
    max_candidates          = cfg.Value("DSOFT_params", "max_candidates");
    num_nz_bins             = cfg.Value("DSOFT_params", "num_nz_bins");
    ignore_lower            = cfg.Value("DSOFT_params", "ignore_lower");

    // GACT first tile
    first_tile_size            = cfg.Value("GACT_first_tile", "first_tile_size");
    first_tile_score_threshold = cfg.Value("GACT_first_tile", "first_tile_score_threshold");

    // GACT extend
    tile_size    = cfg.Value("GACT_extend", "tile_size");
    tile_overlap = cfg.Value("GACT_extend", "tile_overlap");

    // Multi-threading
    num_threads = cfg.Value("Multithreading", "num_threads");

    seed_shape = seed_shape_str.c_str();
    // Ignore lower case in ntcoding
    if (ignore_lower) {
        SetIgnoreLower();
    }

    initRCMap();

    int shape_size = seed_shape.length();

    std::string reference_filename(argv[1]);
    std::string reads_filename(argv[2]);
#ifdef BATCH
    NUM_BLOCKS = std::stoi(argv[3], nullptr);
    THREADS_PER_BLOCK = std::stoi(argv[4], nullptr);
    BATCH_SIZE = NUM_BLOCKS * THREADS_PER_BLOCK;
#endif

#ifdef BATCH
#ifdef GPU
    printf("Using BATCH GPU, batch_size: %d * %d = %d", NUM_BLOCKS, THREADS_PER_BLOCK, BATCH_SIZE);
#else
    printf("Using BATCH, batch_size: %d * %d = %d", NUM_BLOCKS, THREADS_PER_BLOCK, BATCH_SIZE);
#endif	// end GPU
#else
    std::cout << "Running on cpu" << std::endl;
#endif	// end BATCH

    int num_kmer = num_seeds;
    int kmer_count_threshold = dsoft_threshold;

    // LOAD REFERENCE
    std::cout << "\nLoading reference genome ...\n";
    gettimeofday(&start, NULL);

    ParseFastaFile(reference_filename, reference_descrips, reference_seqs, reference_lengths, reference_fileposs);

    reference_string = "";

    int curr_bin = 0;

    for (int i=0; i < reference_seqs.size(); i++) {
        chr_id_to_start_bin[i] =  curr_bin;
        reference_string += reference_seqs[i];
        for (int j = 0; j < (reference_seqs[i].length() / bin_size); j++) {
            bin_to_chr_id[curr_bin++] = i;
        }
        if (reference_seqs[i].length() % bin_size > 0) {
            reference_string += std::string((bin_size - (reference_seqs[i].length() % bin_size)), 'N');
            bin_to_chr_id[curr_bin++] = i;
        }
        std::cout << reference_descrips[i][0] << " length: " << reference_lengths[i] << std::endl;
    }

    reference_length = reference_string.length();
    std::cout << "Reference length: " << (unsigned int) reference_length << std::endl;


    gettimeofday(&end_time, NULL);
    long useconds = end_time.tv_usec - start.tv_usec;
    long seconds = end_time.tv_sec - start.tv_sec;
    long mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time elapsed (loading reference genome): " << mseconds <<" msec" << std::endl;

    // LOAD READS
    std::cout << "\nLoading reads ...\n";
    gettimeofday(&start, NULL);

    ParseFastaFile(reads_filename, reads_descrips, reads_seqs, reads_lengths, reads_fileposs);

    int num_reads = reads_seqs.size();
    for (int i = 0; i < num_reads; i++) {
        std::string rev_read = RevComp(reads_seqs[i]);
        rev_reads_seqs.push_back(rev_read);
    }

    std::cout << "Number of reads: " << num_reads << std::endl;

    gettimeofday(&end_time, NULL);
    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time elapsed (loading reads): " << mseconds <<" msec" << std::endl;

    char* reference_char = (char*) reference_string.c_str();

    reads_char = new char*[num_reads];
    rev_reads_char = new char*[num_reads];
    for (int i =0; i < num_reads; i++) {
        reads_char[i] = (char*) reads_seqs[i].c_str();
        rev_reads_char[i] = (char*) rev_reads_seqs[i].c_str();
    }

    // CONSTRUCT SEED POSITION TABLE
    std::cout << "\nConstructing seed position table ...\n";
    gettimeofday(&start, NULL);

    sa = new SeedPosTable(reference_char, reference_length, seed_shape, seed_occurence_multiple, bin_size);

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time elapsed (seed position table construction): " << mseconds <<" msec" << std::endl;

    uint64_t* seed_offset_vector = new uint64_t[num_kmer];
    uint32_t index = 0;
    int last_N_pos = -1;
    uint64_t sum = 0;
    uint64_t total_num_seeds = 0;
    char nt;

    // RUN D-SOFT TO MAP READS
    gettimeofday(&start, NULL);

    std::cout << "\nFinding candidate bin locations for each read: " << std::endl;

    std::vector<std::thread> align_threads;
    int reads_per_thread = ceil(1.0*num_reads/num_threads);
    for (int k = 0; k < num_reads; k+=reads_per_thread) {
        int last_read = (k+reads_per_thread > num_reads) ? num_reads : k+reads_per_thread;
        align_threads.push_back(std::thread(AlignReads, k, last_read));
    }
    std::cout << align_threads.size() << " threads created\n";
    std::cout << "Synchronizing all threads...\n";
    for (auto& th : align_threads) th.join();

    gettimeofday(&end_time, NULL);

    useconds = end_time.tv_usec - start.tv_usec;
    seconds = end_time.tv_sec - start.tv_sec;
    mseconds = ((seconds) * 1000 + useconds/1000.0) + 0.5;
    std::cout << "Time elapsed (seed table querying): " << mseconds <<" msec" << std::endl;

    return 0;
}

