

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <queue>
#include "align.h"

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

extern int gap_open;
extern int gap_extend;

void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold) {

    std::queue<int> BT_states;

    //output of the function
    std::string aligned_ref_str = "";
    std::string aligned_query_str = "";
    
    int ref_tile_length = tile_size;
    int query_tile_length = tile_size;
    
    int abpos, bbpos;

    static int rpair = 0;

    // beginning for the left or the right extension
    //original starting point for the whole query
    int rev_ref_pos = ref_pos;
    int rev_query_pos = query_pos;
    
    int i = 0;
    int j = 0;
    
    int first_tile_score = 0;
    bool first_tile = true;
   
    printf("GACT: ref_pos: %d, query_pos: %d, ref_length: %d, query_length: %d\n", ref_pos, query_pos, ref_length, query_length);


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
        //printf("after tb: rpair: %d, ref_pos: %d, query_pos: %d, i: %d, j: %d\n", rpair, ref_pos, query_pos, i, j);
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
        //printf("end of while iteration, ref_pos: %d, query_pos: %d, j: %d, i: %d\n", ref_pos, query_pos, j, i);
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
            total_score += sub_mat[5*NtChar2Int(query_nt) + NtChar2Int(ref_nt)];
            open = true;
        }
    }

    printf("End of GACT, rpair: %d, ab: %d, ae: %d, bb: %d, be: %d, score: %d\n", rpair, abpos, ref_pos, bbpos, query_pos, total_score);
    rpair++;

    //std::cout << aligned_ref_str << std::endl << aligned_query_str << std::endl;
    //std::cout << "First tile score: " <<  first_tile_score << " Total score: " << total_score << std::endl << std::endl;
} // end GACT()

void GACT_Batch(std::vector<std::string> ref_strs, \
    std::vector<std::string> query_strs, \
    int tile_size, int tile_overlap, \
    int first_tile_score_threshold, \
    int ref_pos, int query_pos, int num_blocks, int threads_per_block){

    int num_rpairs = ref_strs.size();
    int early_terminate = tile_size - tile_overlap;

    std::vector<int> ref_poss(num_rpairs, ref_pos);
    std::vector<int> query_poss(num_rpairs, query_pos);

    /*for(int i = 0; i < num_rpairs; ++i){
        ref_poss.push_back(ref_pos);
        query_poss.push_back(query_pos);
    }*/

    std::vector<std::queue<int> > BT_statess(num_rpairs);
   
    //output of the function
    std::vector<std::string> aligned_ref_strs(num_rpairs);
    std::vector<std::string> aligned_query_strs(num_rpairs);
    
    //length of the complete sequences
    std::vector<int> ref_lengths(num_rpairs);
    std::vector<int> query_lengths(num_rpairs);
    for(int i = 0; i < num_rpairs; ++i){
        ref_lengths[i] = ref_strs[i].length();
        query_lengths[i] = query_strs[i].length();
    }
    
    std::vector<int> abpos(num_rpairs);
    std::vector<int> bbpos(num_rpairs);

    // beginning for the left or the right extension
    //original starting point for the whole query
    std::vector<int> rev_ref_poss = ref_poss;
    std::vector<int> rev_query_poss = query_poss;

    std::vector<int> first_tile_scores(num_rpairs);

    std::vector<char> rpair_statuss(num_rpairs, 0);     // 0: reverse, 1: forward
    std::vector<char> firsts(num_rpairs, 1);            // 1: first tile

    int BATCH_SIZE = num_blocks * threads_per_block;
    int rpairidx = BATCH_SIZE;
    int match = 1;
    int mismatch = -1;
    int rpairsdone = 0;

    std::vector<int> assignments(BATCH_SIZE);
    std::vector<int> forward(BATCH_SIZE);    
    std::vector<int> ref_tile_lengths(BATCH_SIZE);
    std::vector<int> query_tile_lengths(BATCH_SIZE);
    if(num_rpairs < BATCH_SIZE){printf("ERROR not enough rpairs for BATCH\n");return;}
    for(int i = 0; i < BATCH_SIZE; ++i){
        assignments[i] = i;
    }

    std::vector<std::string> ref_seqs(BATCH_SIZE);
    std::vector<std::string> query_seqs(BATCH_SIZE);
    std::vector<int> ref_lens(BATCH_SIZE);
    std::vector<int> query_lens(BATCH_SIZE);
    //std::vector<int> ref_poss_b(BATCH_SIZE);
    //std::vector<int> query_poss_b(BATCH_SIZE);
    std::vector<char> terminate(BATCH_SIZE);
    std::vector<char> reverses(BATCH_SIZE);
    std::vector<char> firsts_b(BATCH_SIZE);

#ifdef TIME
    std::chrono::high_resolution_clock::time_point t1, t2;
    auto time_gpu = 0, time_loop = 0;
    t1 = std::chrono::high_resolution_clock::now();
#endif

#ifdef GPU
    GPU_storage s;
    GPU_init(BATCH_SIZE, tile_size, tile_overlap, gap_open, gap_extend, match, mismatch, early_terminate, &s);
#endif
int batch_no = 0;
    while(rpairsdone < num_rpairs){
//printf("batch_no: %d\n", batch_no++);
        for(int j = 0; j < BATCH_SIZE; ++j){
            int next_forward = 0;
            int next_rpair = 0;
            int rpair = assignments[j];

            if(rpair == -1){continue;}

            int ref_pos = ref_poss[rpair];
            int query_pos = query_poss[rpair];
            int ref_length = ref_lengths[rpair];
            int query_length = query_lengths[rpair];

            // prepare assignments
            if(rpair_statuss[rpair] == 0){
                if(ref_pos <= 0 || query_pos <= 0 || terminate[j]){
                    //printf("T%d reverse dir done\n", j);
                    abpos[rpair] = ref_pos;
                    bbpos[rpair] = query_pos;
                    ref_pos = rev_ref_poss[rpair];
                    query_pos = rev_query_poss[rpair];
                    ref_poss[rpair] = ref_pos;
                    query_poss[rpair] = query_pos;
                    ref_poss[j] = ref_pos;
                    query_poss[j] = query_pos;
                    next_forward = 1;
                    rpair_statuss[rpair] = 1;
                    terminate[j] = 0;
                }
            }else{
                if(ref_pos >= ref_length || query_pos >= query_length || terminate[j]){
                    printf("End of GACT, rpair: %d, ab: %d, bb: %d, ae: %d, be: %d\n", rpair, abpos[rpair], bbpos[rpair], ref_pos, query_pos);
                    next_rpair = 1;
                    rpairsdone++;
                    assignments[j] = rpairidx;
                    if(rpairidx >= num_rpairs){
                        assignments[j] = -1;
                        //printf("T%d idle\n", j);
                        continue;
                    }
                    rpair = rpairidx++;
                    terminate[j] = 0;
                    ref_pos = ref_poss[rpair];
                    query_pos = query_poss[rpair];
                }
            }

            // prepare batch
            // if first tile
            firsts_b[j] = firsts[rpair];

            // if reverse
            if(rpair_statuss[rpair] == 0){
                ref_lens[j] = (ref_pos > tile_size) ? tile_size : ref_pos;
                query_lens[j] = (query_pos > tile_size) ? tile_size : query_pos;
                ref_seqs[j] = ref_strs[rpair].substr(ref_pos-ref_lens[j], ref_lens[j]);
                query_seqs[j] = query_strs[rpair].substr(query_pos-query_lens[j], query_lens[j]);
                reverses[j] = 0;
            }else{      // else forward
                ref_lens[j] = (ref_pos + tile_size < ref_length) ? tile_size : ref_length - ref_pos;
                query_lens[j] = (query_pos + tile_size < query_length) ? tile_size : query_length - query_pos;
                ref_seqs[j] = ref_strs[rpair].substr(ref_pos, ref_lens[j]);
                query_seqs[j] = query_strs[rpair].substr(query_pos, query_lens[j]);
                reverses[j] = 1;
            }
            //printf("T%d assignment rpair: %d, reverse: %d, first: %d, ref_pos: %d, query_pos: %d\n", j, rpair, reverses[j], firsts_b[j], ref_pos, query_pos);
        }   // end prepare batch

#ifdef TIME
        t2 = std::chrono::high_resolution_clock::now();
        time_loop += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        t1 = std::chrono::high_resolution_clock::now();
#endif

#ifdef GPU
        BT_statess = Align_Batch_GPU(ref_seqs, query_seqs, ref_lens, query_lens, sub_mat, gap_open, gap_extend, ref_lens, query_lens, reverses, firsts_b, early_terminate, tile_size, &s, num_blocks, threads_per_block);
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
            int rpair = assignments[t];

            if(rpair == -1){continue;}

            int i = 0;
            int j = 0;

            std::queue<int> BT_states = BT_statess[t];
            bool first_tile = firsts[rpair];
            int ref_pos = ref_poss[rpair];
            int query_pos = query_poss[rpair];
            int ref_tile_length = ref_lens[t];
            int query_tile_length = query_lens[t];
            int tile_score = BT_states.front();
            int first_tile_score;
            BT_states.pop();
            
            // if reverse
            if(rpair_statuss[rpair] == 0){
                //printf("T%d tb in reverse dir\n");
                if (first_tile) {
                    ref_pos = ref_pos - ref_tile_length + BT_states.front();
                    BT_states.pop();
                    query_pos = query_pos - query_tile_length + BT_states.front();
                    BT_states.pop();
                    rev_ref_poss[rpair] = ref_pos;
                    rev_query_poss[rpair] = query_pos;
                    first_tile_scores[rpair] = tile_score;
                    if (tile_score < first_tile_score_threshold) {
                        break;
                    }
                }
                while (!BT_states.empty()) {
                    first_tile = false;
                    int state = BT_states.front();
                    BT_states.pop();
                    if (state == M) {
                        aligned_ref_strs[rpair].insert(0, 1, ref_strs[rpair][ref_pos - j - 1]);
                        aligned_query_strs[rpair].insert(0, 1, query_strs[rpair][query_pos - i - 1]);
                        i += 1;
                        j += 1;
                    }
                    if (state == I) {
                        aligned_ref_strs[rpair].insert(0, 1, ref_strs[rpair][ref_pos - j - 1]);
                        aligned_query_strs[rpair].insert(0, 1, '-');
                        j += 1;
                    }
                    if (state == D) {
                        aligned_ref_strs[rpair].insert(0, 1, '-');
                        aligned_query_strs[rpair].insert(0, 1, query_strs[rpair][query_pos - i - 1]);
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
                    first_tile_scores[rpair] = tile_score;
                    if (tile_score < first_tile_score_threshold) {
                        break;
                    }
                }
                while (!BT_states.empty()) {
                    first_tile = false;
                    int state = BT_states.front();
                    BT_states.pop();
                    if (state == M) {
                        aligned_ref_strs[rpair] += ref_strs[rpair][ref_pos + j];
                        aligned_query_strs[rpair] += (query_strs[rpair][query_pos + i]);
                        i += 1;
                        j += 1;
                    }
                    if (state == I) {
                        aligned_ref_strs[rpair] += ref_strs[rpair][ref_pos + j];
                        aligned_query_strs[rpair] += '-';
                        j += 1;
                    }
                    if (state == D) {
                        aligned_ref_strs[rpair] += '-';
                        aligned_query_strs[rpair] += query_strs[rpair][query_pos + i];
                        i += 1;
                    }
                }
            
                ref_pos += (j);
                query_pos += (i);
            }   // end traceback

            firsts[rpair] = first_tile;

            if(i == 0 || j == 0){
                terminate[t] = 1;
            }
            ref_poss[rpair] = ref_pos;
            query_poss[rpair] = query_pos;
            //printf("T%d after tb: rpair: %d, ref_pos: %d, query_pos: %d, terminate: %d, rev: %d, i: %d, j: %d\n", t, rpair, ref_pos, query_pos, terminate[t], rpair_statuss[rpair], i, j);
        } // end postprocess

    } // end main loop

#ifdef TIME
        t2 = std::chrono::high_resolution_clock::now();
        time_loop += std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        printf("time_loop: %d ms, time_gpu: %d ms\n", time_loop, time_gpu);
#endif

} // end GACT_Batch()



