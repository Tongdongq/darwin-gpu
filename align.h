#include <iostream>
#include <queue>
#include "ntcoding.h"

#define NOCUDA
#include "cuda_header.h"
#undef NOCUDA

#define INF (1 << 30)
#define MAX_TILE_SIZE 2049

#define MATCH 1
#define MISMATCH -1
#define GAP_OPEN -1
#define GAP_EXTEND -1

typedef int AlnOp;
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};


#ifdef GPU
void GPU_init(int BATCH_SIZE, int tile_size, int tile_overlap, int gap_open, int gap_extend, int match, int mismatch, int early_terminate, GPU_storage *s);
std::vector<std::queue<int> > Align_Batch_GPU(std::vector<std::string> ref_seqs, std::vector<std::string> query_seqs, std::vector<int> ref_lens, std::vector<int> query_lens, int *sub_mat, int gap_open, int gap_extend, std::vector<int> ref_poss, std::vector<int> query_poss, std::vector<char> reverses, std::vector<char> firsts, int early_terminate, int tile_size, GPU_storage *s, int num_blocks, int threads_per_block);
#endif

using namespace std;

std::queue<int> AlignWithBT(char* ref_seq, long long int ref_len,char* query_seq, long long int query_len, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate);

std::vector<std::queue<int> > Align_Batch(std::vector<std::string> ref_seqs, std::vector<std::string> query_seqs, std::vector<int> ref_lens, std::vector<int> query_lens, int *sub_mat, int gap_open, int gap_extend, std::vector<int> ref_poss_b, std::vector<int> query_poss_b, std::vector<char> reverses, std::vector<char> firsts, int early_terminate);
