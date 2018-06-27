
#ifndef GACT_H
#define GACT_H

#ifdef GPU
    #define BATCH 1
#endif

// actually declared in reference_guided.cpp, but gact.cpp also needs them
#ifdef BATCH
    extern int NUM_BLOCKS;
    extern int THREADS_PER_BLOCK;
    extern int BATCH_SIZE;
#endif
extern int tile_size;
extern int tile_overlap;
extern int first_tile_score_threshold;


typedef struct {
    int ref_id;     // id of reference read
    int query_id;   // id of query read
    int ref_pos;    // start of next tile
    int query_pos;
    int ref_bpos;   // indicates where forward dir should start (during reverse), indicates begin of alignment (during forward)
    int query_bpos;
    int score;      // current score for current GACT_call
    int first_tile_score;     // score of first tile
    char first;     // 1: first tile of the GACT_call, 0: not first tile
    char reverse;   // 1: operate in reverse (towards pos = 0), 0: operate in forward direction (towards pos = len)
} GACT_call;

// each GPU thread in the batch has its own GACT_call assignment
typedef struct {
    int callidx;
} GPU_assign;

// implemented in gact.cpp
void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold);

// implemented in gact.cpp
void GACT_Batch(GACT_call *calls, int num_calls, bool complement);

#endif  // GACT_H
