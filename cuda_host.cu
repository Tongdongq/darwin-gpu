
#include <vector>
#include <queue>
#include <string>

#ifdef GPU

#include "cuda_header.h"
#include "gact.h"

std::vector<std::queue<int> > Align_Batch_GPU(std::vector<std::string> ref_seqs, std::vector<std::string> query_seqs, std::vector<int> ref_lens, std::vector<int> query_lens, int *sub_mat, int gap_open, int gap_extend, std::vector<int> ref_poss, std::vector<int> query_poss, std::vector<char> reverses, std::vector<char> firsts, int early_terminate, int tile_size, GPU_storage *s, int num_blocks, int threads_per_block){

  std::vector<std::queue<int> > result;

  int BATCH_SIZE = num_blocks * threads_per_block;

 /* for(int j = 0; j < BATCH_SIZE; ++j){
    char *ref_seq = (char*)ref_seqs[j].c_str();
    char *query_seq = (char*)query_seqs[j].c_str();
    long long int ref_len = ref_lens[j];
    long long int query_len = query_lens[j];
    int ref_pos = ref_poss_b[j];
    int query_pos = query_poss_b[j];
    bool reverse = (reverses[j] == 1);
    bool first = (firsts[j] == 1);

    std::queue<int> BT_states;

    BT_states = AlignWithBT(ref_seq, ref_len, query_seq, query_len, \
      sub_mat, gap_open, gap_extend, query_pos, ref_pos, reverse, first, early_terminate);

    result.push_back(BT_states);
  }//*/

  char *ref_seqs_b;
  char *query_seqs_b;
  int *ref_lens_b;
  int *query_lens_b;
  int *ref_poss_b;
  int *query_poss_b;
  char *reverses_b;
  char *firsts_b;
  const char *ref_seqs_d = s->ref_seqs_d;
  const char *query_seqs_d = s->query_seqs_d;
  const int *ref_lens_d = s->ref_lens_d;
  const int *query_lens_d = s->query_lens_d;
  const int *ref_poss_d = s->ref_poss_d;
  const int *query_poss_d = s->query_poss_d;
  const char *reverses_d = s->reverses_d;
  const char *firsts_d = s->firsts_d;
  int *outs_b, *outs_d = s->outs_d;

  int ref_curr = 0, query_curr = 0;

  ref_seqs_b = (char*)malloc(BATCH_SIZE * tile_size);
  query_seqs_b = (char*)malloc(BATCH_SIZE * tile_size);
  ref_lens_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  query_lens_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  ref_poss_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  query_poss_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  reverses_b = (char*)malloc(BATCH_SIZE);
  firsts_b = (char*)malloc(BATCH_SIZE);
  outs_b = (int*)malloc(BATCH_SIZE * sizeof(int) * 2 * tile_size);

  for(int t = 0; t < BATCH_SIZE; ++t){
#ifdef COALESCE_BASES
    for(int j = 0; j < ref_lens[t]; ++j){
      ref_seqs_b[t+j*BATCH_SIZE] = ref_seqs[t].c_str()[j];
    }
    for(int j = 0; j < query_lens[t]; ++j){
      query_seqs_b[t+j*BATCH_SIZE] = query_seqs[t].c_str()[j];
    }
#else
    memcpy(ref_seqs_b + ref_curr, ref_seqs[t].c_str(), ref_lens[t]);
    memcpy(query_seqs_b + query_curr, query_seqs[t].c_str(), query_lens[t]);
#endif
    ref_curr += tile_size;
    query_curr += tile_size;
    ref_lens_b[t] = ref_lens[t];
    query_lens_b[t] = query_lens[t];
    ref_poss_b[t] = ref_poss[t];
    query_poss_b[t] = query_poss[t];
    reverses_b[t] = reverses[t];
    firsts_b[t] = firsts[t];
  }



  cudaSafeCall(cudaMemcpy((void*)ref_seqs_d, ref_seqs_b, ref_curr, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_seqs_d, query_seqs_b, query_curr, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)ref_lens_d, ref_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_lens_d, query_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)ref_poss_d, ref_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_poss_d, query_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)reverses_d, reverses_b, BATCH_SIZE, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)firsts_d, firsts_b, BATCH_SIZE, cudaMemcpyHostToDevice));

  int NUM_BLOCKS = num_blocks;
  int THREADSPERBLOCK = threads_per_block;

  Align_Kernel<<<NUM_BLOCKS, THREADSPERBLOCK>>>(ref_seqs_d, query_seqs_d, \
    ref_lens_d, query_lens_d, ref_poss_d, query_poss_d, \
    reverses_d, firsts_d, outs_d, s->matricess_d);

  cudaSafeCall(cudaDeviceSynchronize());

  //printf("kernel done\n");

  cudaSafeCall(cudaMemcpy(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 2 * tile_size, cudaMemcpyDeviceToHost));

  std::queue<int> BT_states;
  int off = 0;

  for(int t = 0; t < BATCH_SIZE; ++t){
    BT_states = std::queue<int>();
    for(int i = 1; i <= outs_b[off]; ++i){
      BT_states.push(outs_b[i+off]);
    }
    off += (2 * tile_size);
    result.push_back(BT_states);
  }

  return result;
}


void GPU_init(int tile_size, int tile_overlap, int gap_open, int gap_extend, int match, int mismatch, int early_terminate, GPU_storage *s)
{
  cudaSafeCall(cudaSetDevice(0));         // select Tesla K40c on CE-cuda server

  cudaSafeCall(cudaMemcpyToSymbol(_tile_size, &(tile_size), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_tile_overlap, &(tile_overlap), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_gap_open, &(gap_open), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_gap_extend, &(gap_extend), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_match, &(match), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_mismatch, &(mismatch), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_early_terminate, &(early_terminate), sizeof(int), 0, cudaMemcpyHostToDevice));

  cudaSafeCall(cudaMalloc((void**)&(s->ref_seqs_d), BATCH_SIZE*tile_size));
  cudaSafeCall(cudaMalloc((void**)&(s->query_seqs_d), BATCH_SIZE*tile_size));
  cudaSafeCall(cudaMalloc((void**)&(s->ref_lens_d), BATCH_SIZE*sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&(s->query_lens_d), BATCH_SIZE*sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&(s->ref_poss_d), BATCH_SIZE*sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&(s->query_poss_d), BATCH_SIZE*sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&(s->reverses_d), BATCH_SIZE));
  cudaSafeCall(cudaMalloc((void**)&(s->firsts_d), BATCH_SIZE));
  cudaSafeCall(cudaMalloc((void**)&(s->outs_d), BATCH_SIZE*sizeof(int)*2*tile_size));
  cudaSafeCall(cudaMalloc((void**)&(s->matricess_d), BATCH_SIZE*sizeof(int)*(tile_size+1)*(tile_size+9)));

  // set print buffer size (debug only)
  printf("NOTE: print buffer size is made larger\n");
  cudaSafeCall(cudaDeviceSetLimit(cudaLimitPrintfFifoSize, (1<<20)*50));
}

#endif // GPU
