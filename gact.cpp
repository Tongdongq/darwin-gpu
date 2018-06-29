

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <queue>
#include <mutex>
#include "align.h"
#include "gact.h"

#ifdef TIME
    #include <chrono>
#endif

int sub_mat[25] = {
    1, -1, -1, -1, 0,
    -1, 1, -1, -1, 0,
    -1, -1, 1, -1, 0,
    -1, -1, -1, 1, 0,
    0, 0, 0, 0, 0
};

int match_score = 1;
int mismatch_score = -1;

// declared in cuda_header.h
extern int gap_open;
extern int gap_extend;

// declared in reference_guided.cpp
extern std::vector<std::string> reference_seqs;
extern std::vector<long long int> reference_lengths;
extern std::vector<std::string> reads_seqs;
extern std::vector<std::string> rev_reads_seqs;
extern std::vector<long long int> reads_lengths;

std::mutex io_lock;

void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold, int ref_id, int query_id) {

    std::queue<int> BT_states;

    //output of the function
    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    
    int abpos, bbpos;

    static int callidx = 0;

    // beginning for the left or the right extension
    //original starting point for the whole query
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
    
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;
   
    //printf("GACT: ref_pos: %d, query_pos: %d, ref_length: %d, query_length: %d\n", ref_pos, query_pos, ref_length, query_length);


    // not for the first tile
    while ((ref_pos > 0) && (query_pos > 0) && (((i > 0) && (j > 0)) || first_tile)) {
        //change the tile length if elements less than that of the tile size
        ref_tile_length = (ref_pos > tile_size) ? tile_size : ref_pos;
        query_tile_length = (query_pos > tile_size) ? tile_size : query_pos;

        BT_states = AlignWithBT(\
            ref_str + ref_pos - ref_tile_length, \
            ref_tile_length, \
            query_str + query_pos - query_tile_length, \
            query_tile_length, \
            sub_mat, gap_open, gap_extend, \
            query_tile_length, ref_tile_length, false, \
            first_tile, (tile_size - tile_overlap));
        i = 0;
        j = 0;
        int tile_score = BT_states.front();
        BT_states.pop();

        if (first_tile) {
            ref_pos = ref_pos - ref_tile_length + BT_states.front();
            BT_states.pop();
            query_pos = query_pos - query_tile_length + BT_states.front();
            BT_states.pop();
            rev_ref_pos = ref_pos;
            rev_query_pos = query_pos;
            first_tile_score = tile_score;
            if (tile_score < first_tile_score_threshold) {
                break;
            }
        }
        while (!BT_states.empty()) {
            first_tile = false;
            int state = BT_states.front();
            BT_states.pop();
            if (state == M) {
                aligned_ref_str.insert(0, 1, ref_str[ref_pos - j - 1]);
                aligned_query_str.insert(0, 1, query_str[query_pos - i - 1]);
                i += 1;
                j += 1;
            }
            if (state == I) {
                aligned_ref_str.insert(0, 1, ref_str[ref_pos - j - 1]);
                aligned_query_str.insert(0, 1, '-');
                j += 1;
            }
            if (state == D) {
                aligned_ref_str.insert(0, 1, '-');
                aligned_query_str.insert(0, 1, query_str[query_pos - i - 1]);
                i += 1;
            }
        }
        ref_pos -= (j);
        query_pos -= (i);
        //printf("after tb: callidx: %d, ref_pos: %d, query_pos: %d, i: %d, j: %d, rev: 1\n", callidx, ref_pos, query_pos, i, j);
    }

    abpos = ref_pos;
    bbpos = query_pos;
    ref_pos = rev_ref_pos;
    query_pos = rev_query_pos;
    //printf("rev wave done, start with for wave\n");
    i = tile_size;
    j = tile_size;
    
    //starts with the first tile
    while ((ref_pos < ref_length) && (query_pos < query_length) && (((i > 0) && (j > 0)) || first_tile)) {
        
        ref_tile_length = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
        query_tile_length = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
        
        BT_states = AlignWithBT (\
            ref_str + ref_pos, \
            ref_tile_length, \
            query_str + query_pos, \
            query_tile_length, \
            sub_mat, gap_open, gap_extend, \
            query_tile_length, ref_tile_length, true, \
            first_tile, (tile_size - tile_overlap));

        i = 0;
        j = 0;
        int tile_score = BT_states.front();
        BT_states.pop();

        if (first_tile) {
            ref_pos = ref_pos + ref_tile_length - BT_states.front();
            BT_states.pop();
            query_pos = query_pos + query_tile_length - BT_states.front();
            BT_states.pop();
            first_tile_score = tile_score;
            if (tile_score < first_tile_score_threshold) {
                break;
            }
        }
        while (!BT_states.empty()) {
            first_tile = false;
            int state = BT_states.front();
            BT_states.pop();
            if (state == M) {
                aligned_ref_str += ref_str[ref_pos + j];
                aligned_query_str += (query_str[query_pos + i]);
                i += 1;
                j += 1;
            }
            if (state == I) {
                aligned_ref_str += ref_str[ref_pos + j];
                aligned_query_str += '-';
                j += 1;
            }
            if (state == D) {
                aligned_ref_str += '-';
                aligned_query_str += query_str[query_pos + i];
                i += 1;
            }
        }
    
        ref_pos += (j);
        query_pos += (i);
        //printf("after tb: callidx: %d, ref_pos: %d, query_pos: %d, i: %d, j: %d, rev: 0\n", callidx, ref_pos, query_pos, i, j);
    }

    int total_score = 0;
    bool open = true;
    for (uint32_t j = 0; j < aligned_ref_str.length(); j++) {
        char ref_nt = aligned_ref_str[j];
        char query_nt = aligned_query_str[j];
        if (ref_nt == '-' || query_nt == '-') {
            total_score += (open) ? gap_open : gap_extend;
            open = false;
        }
        else {
            total_score += (query_nt == ref_nt) ? match_score : mismatch_score;
            open = true;
        }
    }
io_lock.lock();
    //printf("End of GACT, callidx: %d, ab: %d, ae: %d, bb: %d, be: %d, score: %d\n", callidx, abpos, ref_pos, bbpos, query_pos, total_score);
    printf("End of GACT, ref_id: %d, query_id: %d, ab: %d, ae: %d, bb: %d, be: %d, score: %d\n", ref_id, query_id, abpos, ref_pos, bbpos, query_pos, total_score);
io_lock.unlock();
    callidx++;

    //std::cout << aligned_ref_str << std::endl << aligned_query_str << std::endl;
    //std::cout << "First tile score: " <<  first_tile_score << " Total score: " << total_score << std::endl << std::endl;
} // end GACT()

#ifdef BATCH
#ifdef GPU
void GACT_Batch(std::vector<GACT_call> calls, int num_calls, bool complement, int offset, GPU_storage *s)
#else
void GACT_Batch(std::vector<GACT_call> calls, int num_calls, bool complement, int offset)
#endif
{
    int early_terminate = tile_size - tile_overlap;

    std::vector<std::queue<int> > BT_statess(BATCH_SIZE);

    std::vector<std::string> *reads_seqs_p;

    if(complement == true){
        reads_seqs_p = &rev_reads_seqs;
    }else{
        reads_seqs_p = &reads_seqs;
    }

    printf("GACT_Batch, num_calls: %d\n", num_calls);

    //output of the function
    std::vector<std::string> aligned_ref_strs(num_calls);
    std::vector<std::string> aligned_query_strs(num_calls);

    std::vector<int> first_tile_scores(num_calls);

    int next_callidx = BATCH_SIZE;
    //int next_callidx = 5000;
    int calls_done = 0;

    std::vector<int> assignments(BATCH_SIZE);
    std::vector<int> ref_tile_lengths(BATCH_SIZE);
    std::vector<int> query_tile_lengths(BATCH_SIZE);
    if(num_calls < BATCH_SIZE){
        printf("WARNING not enough callidxs for BATCH\n");
        int i = 0;
        for(; i < num_calls; ++i){
            assignments[i] = i;
        }
        for(; i < BATCH_SIZE; ++i){
            assignments[i] = -1;
        }
    }else{
        for(int i = 0; i < BATCH_SIZE; ++i){
            assignments[i] = i;
        }
    }

    // terminate: 0: continue, 1: stop due to low number of tb steps, 2: stop due to low first tile score
    std::vector<char> terminate(BATCH_SIZE);
    std::vector<std::string> ref_seqs(BATCH_SIZE);
    std::vector<std::string> query_seqs(BATCH_SIZE);
    std::vector<int> ref_lens(BATCH_SIZE);
    std::vector<int> query_lens(BATCH_SIZE);
    std::vector<char> reverses(BATCH_SIZE);
    std::vector<char> firsts_b(BATCH_SIZE);

#ifdef TIME
    std::chrono::high_resolution_clock::time_point t1, t2;
    auto time_gpu = 0, time_loop = 0;
    t1 = std::chrono::high_resolution_clock::now();
#endif

int batch_no = 0;
    while(calls_done < num_calls){
//printf("batch_no: %d\n", batch_no++);
        for(int t = 0; t < BATCH_SIZE; ++t){
            char next_call = 0;
            int callidx = assignments[t];
            GACT_call *c = &(calls[callidx]);

            if(callidx == -1){
                ref_lens[t] = -1;
                continue;
            }

            int ref_pos = c->ref_pos;
            int query_pos = c->query_pos;
            int ref_length = reference_lengths[c->ref_id];
            int query_length = reads_lengths[c->query_id];

            // prepare assignments
            if(c->reverse == 1){
                if(ref_pos <= 0 || query_pos <= 0 || terminate[t]){
                    if(terminate[t] == 2){
                        next_call = 1;
                    }else{
                        //printf("T%d reverse dir done\n", t);
                        // store begin of alignment in ref_bpos and query_bpos
                        //printf("T%d reverse done, ref_pos: %d, query_pos: %d\n", t, ref_pos, query_pos);
                        int t1 = c->ref_bpos;
                        int t2 = c->query_bpos;
                        c->ref_bpos = ref_pos;
                        c->query_bpos = query_pos;
                        ref_pos = t1;
                        query_pos = t2;
                        c->ref_pos = t1;
                        c->query_pos = t2;
                        c->reverse = 0;
                        terminate[t] = 0;
                    }
                }
            }else{
                if(ref_pos >= ref_length || query_pos >= query_length || terminate[t]){
                    if(terminate[t] != 2){
io_lock.lock();
                        printf("End of GACT, ref_id: %d, query_id: %d, ab: %d, ae: %d, bb: %d, be: %d\n", c->ref_id, c->query_id, c->ref_bpos, ref_pos, c->query_bpos, query_pos);
io_lock.unlock();
                        calls_done++;
                        assignments[t] = next_callidx;
                        if(next_callidx >= num_calls){
                            assignments[t] = -1;
                            ref_lens[t] = -1;
                            //printf("T%d idle\n", t);
                            continue;
                        }
                    }
                    next_call = 1;
                }
            }

            if(next_call == 1){                
                callidx = next_callidx++;
                c = &(calls[callidx]);
                ref_pos = c->ref_pos;
                query_pos = c->query_pos;
                ref_length = reference_lengths[c->ref_id];
                query_length = reads_lengths[c->query_id];
                terminate[t] = 0;
            }

            // prepare batch
            // if first tile
            firsts_b[t] = c->first;

            // if reverse
            if(c->reverse == 1){
                ref_lens[t] = (ref_pos > tile_size) ? tile_size : ref_pos;
                query_lens[t] = (query_pos > tile_size) ? tile_size : query_pos;
                ref_seqs[t] = reference_seqs[c->ref_id].substr(ref_pos-ref_lens[t], ref_lens[t]);
                query_seqs[t] = reads_seqs_p->at(c->query_id).substr(query_pos-query_lens[t], query_lens[t]);
                reverses[t] = 1;
            }else{      // else forward
                ref_lens[t] = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
                query_lens[t] = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
                ref_seqs[t] = reference_seqs[c->ref_id].substr(ref_pos, ref_lens[t]);
                query_seqs[t] = reads_seqs_p->at(c->query_id).substr(query_pos, query_lens[t]);
                reverses[t] = 0;
            }
            //printf("T%d assignment callidx: %d, reverse: %d, first: %d, ref_pos: %d, query_pos: %d\n", t, callidx, reverses[t], firsts_b[t], ref_pos, query_pos);
        }   // end prepare batch

#ifdef TIME
        t2 = std::chrono::high_resolution_clock::now();
        time_loop += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        t1 = std::chrono::high_resolution_clock::now();
#endif

#ifdef GPU
        BT_statess = Align_Batch_GPU(ref_seqs, query_seqs, ref_lens, query_lens, sub_mat, gap_open, gap_extend, ref_lens, query_lens, reverses, firsts_b, early_terminate, tile_size, s, NUM_BLOCKS, THREADS_PER_BLOCK);
#else
        //BT_statess = Align_Batch(ref_seqs, query_seqs, ref_lens, query_lens, sub_mat, gap_open, gap_extend, ref_poss_b, query_poss_b, reverses, firsts_b, early_terminate);
        BT_statess = Align_Batch(ref_seqs, query_seqs, ref_lens, query_lens, sub_mat, gap_open, gap_extend, ref_lens, query_lens, reverses, firsts_b, early_terminate);
#endif

#ifdef TIME
        t2 = std::chrono::high_resolution_clock::now();
        time_gpu += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        t1 = std::chrono::high_resolution_clock::now();
#endif
        // postprocess
        for(int t = 0; t < BATCH_SIZE; ++t){
            int callidx = assignments[t];
            GACT_call *c = &(calls[callidx]);

            if(callidx == -1){continue;}

            int i = 0;
            int j = 0;

            std::queue<int> BT_states = BT_statess[t];
            bool first_tile = c->first;
            int ref_pos = c->ref_pos;
            int query_pos = c->query_pos;
            int ref_tile_length = ref_lens[t];
            int query_tile_length = query_lens[t];
            int tile_score = BT_states.front();
            int first_tile_score;
            BT_states.pop();

            // if reverse
            if(c->reverse == 1){
                //printf("T%d tb in reverse dir\n");
                if (first_tile) {
                    ref_pos = ref_pos - ref_tile_length + BT_states.front();
                    BT_states.pop();
                    query_pos = query_pos - query_tile_length + BT_states.front();
                    BT_states.pop();
                    c->ref_bpos = ref_pos;
                    c->query_bpos = query_pos;
                    c->first_tile_score = tile_score;
                    if (tile_score < first_tile_score_threshold) {
                        //terminate[t] = 2;
                        terminate[t] = 1;
                        c->ref_pos = ref_pos;
                        c->query_pos = query_pos;
                        continue;
                    }
                }
                while (!BT_states.empty()) {
                    first_tile = false;
                    int state = BT_states.front();
                    BT_states.pop();
                    if (state == M) {
                        aligned_ref_strs[callidx].insert(0, 1, reference_seqs[c->ref_id][ref_pos - j - 1]);
                        aligned_query_strs[callidx].insert(0, 1, reads_seqs_p->at(c->query_id)[query_pos - i - 1]);
                        i += 1;
                        j += 1;
                    }
                    if (state == I) {
                        aligned_ref_strs[callidx].insert(0, 1, reference_seqs[c->ref_id][ref_pos - j - 1]);
                        aligned_query_strs[callidx].insert(0, 1, '-');
                        j += 1;
                    }
                    if (state == D) {
                        aligned_ref_strs[callidx].insert(0, 1, '-');
                        aligned_query_strs[callidx].insert(0, 1, reads_seqs_p->at(c->query_id)[query_pos - i - 1]);
                        i += 1;
                    }
                }
                //printf("T%d done with tb in rev dir, i: %d, j: %d\n", t, i, j);
                ref_pos -= (j);
                query_pos -= (i);
            }else{      // else forward
                if (first_tile) {
                    ref_pos = ref_pos + ref_tile_length - BT_states.front();
                    BT_states.pop();
                    query_pos = query_pos + query_tile_length - BT_states.front();
                    BT_states.pop();
                    c->first_tile_score = tile_score;
                    if (tile_score < first_tile_score_threshold) {
                        //terminate[t] = 2;
                        terminate[t] = 1;
                        c->ref_pos = ref_pos;
                        c->query_pos = query_pos;
                        continue;
                    }
                }
                while (!BT_states.empty()) {
                    first_tile = false;
                    int state = BT_states.front();
                    BT_states.pop();
                    if (state == M) {
                        aligned_ref_strs[callidx] += reference_seqs[c->ref_id][ref_pos + j];
                        aligned_query_strs[callidx] += (reads_seqs_p->at(c->query_id)[query_pos + i]);
                        i += 1;
                        j += 1;
                    }
                    if (state == I) {
                        aligned_ref_strs[callidx] += reference_seqs[c->ref_id][ref_pos + j];
                        aligned_query_strs[callidx] += '-';
                        j += 1;
                    }
                    if (state == D) {
                        aligned_ref_strs[callidx] += '-';
                        aligned_query_strs[callidx] += reads_seqs_p->at(c->query_id)[query_pos + i];
                        i += 1;
                    }
                }
            
                ref_pos += (j);
                query_pos += (i);
            }   // end traceback

            c->first = first_tile;

            if(i == 0 || j == 0){
                terminate[t] = 1;
            }
            c->ref_pos = ref_pos;
            c->query_pos = query_pos;
            //printf("T%d after tb: callidx: %d, ref_pos: %d, query_pos: %d, terminate: %d, rev: %d, i: %d, j: %d\n", t, callidx+offset, ref_pos, query_pos, terminate[t], c->reverse, i, j);
        } // end postprocess

    } // end main loop

#ifdef TIME
        t2 = std::chrono::high_resolution_clock::now();
        time_loop += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        printf("time_loop: %d ms, time_gpu: %d ms\n", time_loop, time_gpu);
#endif

} // end GACT_Batch()
#endif  // BATCH


