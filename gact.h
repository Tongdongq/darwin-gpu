

int NUM_BLOCKS;
int THREADS_PER_BLOCK;
int BATCH_SIZE;

typedef struct {
    int ref_id;     // id of reference read
    int query_id;   // id of query read
    int ref_pos;    // start of next tile
    int query_pos;
    int ref_epos;   // indicates end of alignment
    int query_epos;
    int score;      // current score for current GACT_call
    char reverse;   // 1: operate on reverse read (towards pos = 0), 0: operate on forward read (towards pos = len)
    char terminate; // 1: terminate current GACT_call because previous tile did not have enough alignment, 0: continue with next tile
} GACT_call;

// each GPU thread in the batch has its own GACT_call assignment
typedef struct {
  int callidx;
} GPU_assign;

void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold);

