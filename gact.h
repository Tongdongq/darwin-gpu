

int NUM_BLOCKS;
int THREADS_PER_BLOCK;

typedef struct {
	int ref_id;
	int query_id;
	int ref_pos;
	int query_pos;

} GACT_call;

void GACT (char *ref_str, char *query_str, \
    int ref_length, int query_length, \
    int tile_size, int tile_overlap, \
    int ref_pos, int query_pos, int first_tile_score_threshold);

