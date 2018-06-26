

int NUM_BLOCKS;
int THREADS_PER_BLOCK;
int BATCH_SIZE;

typedef struct {
    int ref_id;
    int query_id;
    int ref_pos;
    int query_pos;
} GACT_call;

// each GPU thread in the batch has its own assignment
// current GACT_call idx is stored in rpair_status[GPU_assign[tid].rpair]
typedef struct {
  int rpair;        // current rpair for GPU thread
} GPU_assign;

typedef struct {
    int idx;        // indicates current GACT_call
    int score;      // current score for current GACT_call
    int ref_pos;    // indicates start of tile
    int query_pos;  // indicates start of tile
    int ref_epos;   // indicates end of alignment
    int query_epos; // indicates end of alignment
    char reverse;   // 1: operate on reverse read (towards pos = 0), 0: operate on forward read (towards pos = len)
    char terminate; // 1: terminate current GACT_call because previous tile did not have enough alignment, 0: continue with next tile
} Rpair_status;

// represents a readpair
// is associated with one or more GACT_calls
typedef struct {
  int ref_id;   // id of associated aread
  int query_id; // id of associated bread
  int bidx;     // idx of this rpairs first GACT_call in large GACT_call array
  int eidx;     // idx of this rpairs last GACT_call in large GACT_call array
} Readpair;

void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold);

