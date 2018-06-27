#include <iostream>
#include <queue>
#include "ntcoding.h"


#define INF (1 << 30)
#define MAX_TILE_SIZE 2049

#define MATCH 1
#define MISMATCH -1
#define GAP_OPEN -1
#define GAP_EXTEND -1

typedef int AlnOp;
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};



using namespace std;

std::queue<int> AlignWithBT(char* ref_seq, long long int ref_len,char* query_seq, long long int query_len, int* sub_mat, int gap_open, int gap_extend, int query_pos, int ref_pos, bool reverse, bool first, int early_terminate);

std::vector<std::queue<int> > Align_Batch(std::vector<std::string> ref_seqs, std::vector<std::string> query_seqs, std::vector<int> ref_lens, std::vector<int> query_lens, int *sub_mat, int gap_open, int gap_extend, std::vector<int> ref_poss_b, std::vector<int> query_poss_b, std::vector<char> reverses, std::vector<char> firsts, int early_terminate);
