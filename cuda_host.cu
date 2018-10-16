
/*
Copyright 2018 Yatish Turakhia, Gill Bejerano and William Dally
Copyright 2018 Tong Dong Qiu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <vector>
#include <queue>
#include <string>

#ifdef GPU

#include "cuda_header.h"
#include "gact.h"


int* Align_Batch_GPU( \
    std::vector<std::string> ref_seqs, std::vector<std::string> query_seqs, \
    std::vector<int> ref_lens, std::vector<int> query_lens, \
    int *sub_mat, int gap_open, int gap_extend, \
    std::vector<int> ref_poss, std::vector<int> query_poss, \
    std::vector<char> reverses, std::vector<char> firsts, \
    int early_terminate, int tile_size, \
    GPU_storage *s, int num_blocks, int threads_per_block){

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

  int ref_curr = 0, query_curr = 0;

  ref_seqs_b = (char*)malloc(BATCH_SIZE * (tile_size+1));
  query_seqs_b = (char*)malloc(BATCH_SIZE * (tile_size+1));
  ref_lens_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  query_lens_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  ref_poss_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  query_poss_b = (int*)malloc(BATCH_SIZE * sizeof(int));
  reverses_b = (char*)malloc(BATCH_SIZE);
  firsts_b = (char*)malloc(BATCH_SIZE);
#ifdef GASAL
  int32_t *query_offsets_b = (int32_t*)malloc(BATCH_SIZE * sizeof(int));
  int32_t *ref_offsets_b = (int32_t*)malloc(BATCH_SIZE * sizeof(int));
  int *outs_b, *outs_d = (s->outs_d);
  outs_b = (int*)malloc(BATCH_SIZE * sizeof(int) * 2 * tile_size);
#else
  int *outs_b, *outs_d = s->outs_d;
  outs_b = (int*)malloc(BATCH_SIZE * sizeof(int) * 2 * tile_size);
#endif

  for(int t = 0; t < BATCH_SIZE; ++t){
#ifndef GASAL
#ifdef COALESCE_BASES
    /*for(int j = 0; j < ref_lens[t]; ++j){
      ref_seqs_b[t+j*BATCH_SIZE] = ref_seqs[t].c_str()[j];
    }
    for(int j = 0; j < query_lens[t]; ++j){
      query_seqs_b[t+j*BATCH_SIZE] = query_seqs[t].c_str()[j];
    }//*/
    if(reverses[t] == 1){
      for(int j = 0; j < ref_lens[t]; ++j){
        ref_seqs_b[t+j*BATCH_SIZE] = ref_seqs[t].c_str()[j];
      }
      for(int j = 0; j < query_lens[t]; ++j){
        query_seqs_b[t+j*BATCH_SIZE] = query_seqs[t].c_str()[j];
      }
    }else{
      for(int j = 0; j < ref_lens[t]; ++j){
        ref_seqs_b[t+j*BATCH_SIZE] = ref_seqs[t].c_str()[ref_lens[t]-j-1];
      }
      for(int j = 0; j < query_lens[t]; ++j){
        query_seqs_b[t+j*BATCH_SIZE] = query_seqs[t].c_str()[query_lens[t]-j-1];
      }
    }//*/
    ref_curr += tile_size;
    query_curr += tile_size;
#else
    if(reverses[t] == 1){
        memcpy(ref_seqs_b + ref_curr, ref_seqs[t].c_str(), ref_lens[t]);
        memcpy(query_seqs_b + query_curr, query_seqs[t].c_str(), query_lens[t]);
    }else{
        for(int j = 0; j < ref_lens[t]; ++j){
            ref_seqs_b[ref_curr+j] = ref_seqs[t].c_str()[ref_lens[t]-j-1];
        }
        for(int j = 0; j < query_lens[t]; ++j){
            query_seqs_b[query_curr+j] = query_seqs[t].c_str()[query_lens[t]-j-1];
        }
    }
    ref_curr += tile_size;
    query_curr += tile_size;
#endif // COALESCE_BASES
#endif // GASAL
    ref_lens_b[t] = ref_lens[t];
    //printf("T%d ref_len: %d\n", t, ref_lens_b[t]);
    query_lens_b[t] = query_lens[t];
    ref_poss_b[t] = ref_poss[t];
    query_poss_b[t] = query_poss[t];
    reverses_b[t] = reverses[t];
    firsts_b[t] = firsts[t];
  }

  int NUM_BLOCKS = num_blocks;
  int THREADSPERBLOCK = threads_per_block;


#ifdef STREAM
#ifdef GASAL
  for(int t = 0; t < BATCH_SIZE; ++t){
    if(ref_lens[t] == -1){continue;}
    ref_offsets_b[t] = ref_curr;
    query_offsets_b[t] = query_curr;
#ifndef COALESCE_PACKED_BASES
    if(reverses[t] == 1){
        memcpy(ref_seqs_b + ref_curr, ref_seqs[t].c_str(), ref_lens[t]);
        memcpy(query_seqs_b + query_curr, query_seqs[t].c_str(), query_lens[t]);
    }else{
      //printf("T%d is forward\n", t);
        for(int j = 0; j < ref_lens[t]; ++j){
            ref_seqs_b[ref_curr+j] = ref_seqs[t].c_str()[ref_lens[t]-j-1];
        }
        //printf("T%d query %d \t", t, query_lens[t]);
        for(int j = 0; j < query_lens[t]; ++j){
          //if(j % 8 == 0){printf(" ");}
            //printf("%d", query_seqs[t].c_str()[query_lens[t]-j-1]);
            query_seqs_b[query_curr+j] = query_seqs[t].c_str()[query_lens[t]-j-1];
        }//printf("\n");
    }
    ref_curr += ref_lens[t];
    query_curr += query_lens[t];
    while(ref_curr % 8 != 0){
      ref_seqs_b[ref_curr++] = 4;
    }
    while(query_curr % 8 != 0){
      query_seqs_b[query_curr++] = 5;
    }
#else
    int i, j, remainder;
    if(reverses[t] == 1){
      for(i = 0; i < ref_lens[t] / 8; ++i){
        memcpy(ref_seqs_b + (i * BATCH_SIZE + t) * 8, ref_seqs[t].c_str() + 8*i, 8);
      }
      remainder = ref_lens[t] % 8;
      if(remainder != 0){
        memcpy(ref_seqs_b + (i * BATCH_SIZE + t) * 8, ref_seqs[t].c_str() + 8*i, remainder);
        for(j = remainder; j < 8; ++j){
          ref_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = 4;
        }
      }
      for(i = 0; i < query_lens[t] / 8; ++i){
        memcpy(query_seqs_b + (i * BATCH_SIZE + t) * 8, query_seqs[t].c_str() + 8*i, 8);
      }
      remainder = query_lens[t] % 8;
      if(remainder != 0){
        memcpy(query_seqs_b + (i * BATCH_SIZE + t) * 8, query_seqs[t].c_str() + 8*i, remainder);
        for(j = remainder; j < 8; ++j){
          query_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = 5;
        }
        }
    }else{
        for(i = 0; i < ref_lens[t] / 8; ++i){
          for(j = 0; j < 8; ++j){
            ref_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = ref_seqs[t].c_str()[ref_lens[t]-8*i-j-1];
          }
        }
        remainder = ref_lens[t] % 8;
        if(remainder != 0){
          for(j = 0; j < remainder; ++j){
            ref_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = ref_seqs[t].c_str()[ref_lens[t]-8*i-j-1];;
          }
          for(; j < 8; ++j){
            ref_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = 4;
          }
        }
        for(i = 0; i < query_lens[t] / 8; ++i){
          for(j = 0; j < 8; ++j){
            query_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = query_seqs[t].c_str()[query_lens[t]-8*i-j-1];
          }
        }
        remainder = query_lens[t] % 8;
        if(remainder % 8 != 0){
          for(j = 0; j < remainder; ++j){
            query_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = query_seqs[t].c_str()[query_lens[t]-8*i-j-1];;
          }
          for(; j < 8; ++j){
            query_seqs_b[(i * BATCH_SIZE + t) * 8 + j] = 5;
          }
        }
    }
    ref_curr = BATCH_SIZE * tile_size;
    query_curr = BATCH_SIZE * tile_size;
#endif
  }
  /*for(int i = 0; i < 640; ++i){
    if(i%8==0)printf(" ");
    printf("%d", query_seqs_b[i]);
  }printf("\n");
  for(int i = 0; i < 640; ++i){
    if(i%8==0)printf(" ");
    printf("%d", ref_seqs_b[i]);
  }printf("\n");//*/
  //printf("150\n");
  uint32_t *packed_ref_seqs_d = s->packed_ref_seqs_d;
  uint32_t *packed_query_seqs_d = s->packed_query_seqs_d;
  int32_t *query_offsets_d = s->query_offsets_d;
  int32_t *ref_offsets_d = s->ref_offsets_d;
  int query_batch_tasks_per_thread = (int)ceil((double)query_curr/(8*THREADS_PER_BLOCK*NUM_BLOCKS));
  int target_batch_tasks_per_thread = (int)ceil((double)ref_curr/(8*THREADS_PER_BLOCK*NUM_BLOCKS));
  cudaStream_t stream = s->stream->stream;
  cudaSafeCall(cudaMemcpyAsync((void*)ref_seqs_d, ref_seqs_b, ref_curr, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_seqs_d, query_seqs_b, query_curr, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)ref_lens_d, ref_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_lens_d, query_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)ref_poss_d, ref_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_poss_d, query_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)reverses_d, reverses_b, BATCH_SIZE, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)firsts_d, firsts_b, BATCH_SIZE, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_offsets_d, query_offsets_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)ref_offsets_d, ref_offsets_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));

  gasal_pack_kernel<<<NUM_BLOCKS, THREADSPERBLOCK, 0, stream>>>( \
    (uint32_t*)query_seqs_d, (uint32_t*)ref_seqs_d, \
    packed_query_seqs_d, packed_ref_seqs_d, \
    (uint32_t)query_batch_tasks_per_thread, (uint32_t)target_batch_tasks_per_thread, \
    query_curr / 4, ref_curr / 4);

  cudaSafeCall(cudaStreamSynchronize(stream));

#ifndef GLOBAL
  gasal_local_kernel<<<NUM_BLOCKS, THREADSPERBLOCK, 0, stream>>>( \
    packed_query_seqs_d, packed_ref_seqs_d, \
    query_lens_d, ref_lens_d, \
    query_offsets_d, ref_offsets_d, \
    query_poss_d, ref_poss_d, \
    outs_d, \
    firsts_d, (char*)(s->matrices_d));
#else
  gasal_local_kernel<<<NUM_BLOCKS, THREADSPERBLOCK, 0, stream>>>( \
    packed_query_seqs_d, packed_ref_seqs_d, \
    query_lens_d, ref_lens_d, \
    query_offsets_d, ref_offsets_d, \
    query_poss_d, ref_poss_d, \
    outs_d, \
    firsts_d, (char*)(s->matrices_d), s->global_d);
#endif

  cudaSafeCall(cudaStreamSynchronize(stream));

#ifdef SCORE
  cudaSafeCall(cudaMemcpyAsync(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 2 * tile_size, cudaMemcpyDeviceToHost));
#else
  cudaSafeCall(cudaMemcpyAsync(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 5, cudaMemcpyDeviceToHost));
#endif // SCORE
#else // GASAL
  cudaStream_t stream = s->stream->stream;
  cudaSafeCall(cudaMemcpyAsync((void*)ref_seqs_d, ref_seqs_b, ref_curr, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_seqs_d, query_seqs_b, query_curr, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)ref_lens_d, ref_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_lens_d, query_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)ref_poss_d, ref_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)query_poss_d, query_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)reverses_d, reverses_b, BATCH_SIZE, cudaMemcpyHostToDevice, stream));
  cudaSafeCall(cudaMemcpyAsync((void*)firsts_d, firsts_b, BATCH_SIZE, cudaMemcpyHostToDevice, stream));

  Align_Kernel<<<NUM_BLOCKS, THREADSPERBLOCK, 0, stream>>>(ref_seqs_d, query_seqs_d, \
    ref_lens_d, query_lens_d, ref_poss_d, query_poss_d, \
    reverses_d, firsts_d, outs_d, s->matrices_d);

  cudaSafeCall(cudaStreamSynchronize(stream));

#ifdef SCORE
  cudaSafeCall(cudaMemcpyAsync(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 2 * tile_size, cudaMemcpyDeviceToHost));
#else
  cudaSafeCall(cudaMemcpyAsync(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 5, cudaMemcpyDeviceToHost));
#endif // SCORE
#endif // GASAL
#else // STREAM
  cudaSafeCall(cudaMemcpy((void*)ref_seqs_d, ref_seqs_b, ref_curr, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_seqs_d, query_seqs_b, query_curr, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)ref_lens_d, ref_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_lens_d, query_lens_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)ref_poss_d, ref_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)query_poss_d, query_poss_b, BATCH_SIZE*sizeof(int), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)reverses_d, reverses_b, BATCH_SIZE, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy((void*)firsts_d, firsts_b, BATCH_SIZE, cudaMemcpyHostToDevice));

  Align_Kernel<<<NUM_BLOCKS, THREADSPERBLOCK>>>(ref_seqs_d, query_seqs_d, \
    ref_lens_d, query_lens_d, ref_poss_d, query_poss_d, \
    reverses_d, firsts_d, outs_d, s->matrices_d);

  cudaSafeCall(cudaDeviceSynchronize());

#ifdef SCORE
  cudaSafeCall(cudaMemcpy(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 2 * tile_size, cudaMemcpyDeviceToHost));
#else
  cudaSafeCall(cudaMemcpy(outs_b, outs_d, BATCH_SIZE * sizeof(int) * 5, cudaMemcpyDeviceToHost));
#endif // SCORE
#endif // STREAM


  return outs_b;
}


//void GPU_init(int tile_size, int tile_overlap, int gap_open, int gap_extend, int match, int mismatch, int early_terminate, GPU_storage *s, int num_threads)
void GPU_init(int tile_size, int tile_overlap, int gap_open, int gap_extend, int match, int mismatch, int early_terminate, std::vector<GPU_storage> *s, int num_threads)
{
  cudaSafeCall(cudaSetDevice(0));         // select Tesla K40c on CE-cuda server

  cudaSafeCall(cudaMemcpyToSymbol(_tile_size, &(tile_size), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_tile_overlap, &(tile_overlap), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_gap_open, &(gap_open), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_gap_extend, &(gap_extend), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_match, &(match), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_mismatch, &(mismatch), sizeof(int), 0, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpyToSymbol(_early_terminate, &(early_terminate), sizeof(int), 0, cudaMemcpyHostToDevice));

#ifdef GASAL
  int size_matrices = (tile_size+2)*(tile_size+2);
#ifdef SCORE
  int size_outs = 2*tile_size;
#else
  int size_outs = 5;
#endif // SCORE
#else // GASAL
#ifndef COMPRESS_DIR
  int size_matrices = sizeof(int)*(tile_size+1)*8 + (tile_size+1)*(tile_size+1);
#else
  int size_matrices = sizeof(int)*(tile_size+1)*8 + ((tile_size+1)*(tile_size+1)+1)/2;
#endif // COMPRESS_DIR
#endif // GASAL

  for(int i = 0; i < num_threads; ++i){
    s->push_back(GPU_storage());
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].ref_seqs_d), BATCH_SIZE*(tile_size+1)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].query_seqs_d), BATCH_SIZE*(tile_size+1)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].ref_lens_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].query_lens_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].ref_poss_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].query_poss_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].reverses_d), BATCH_SIZE));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].firsts_d), BATCH_SIZE));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].matrices_d), BATCH_SIZE*size_matrices));
#ifdef GLOBAL
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].global_d), BATCH_SIZE*3*sizeof(short)*tile_size));
#endif
#ifdef STREAM
    (*s)[i].stream = (CUDA_Stream_Holder*)malloc(sizeof(CUDA_Stream_Holder*));
    cudaSafeCall(cudaStreamCreate(&((*s)[i].stream->stream)));
#endif
#ifdef GASAL
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].packed_ref_seqs_d), BATCH_SIZE*tile_size));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].packed_query_seqs_d), BATCH_SIZE*tile_size));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].ref_offsets_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].query_offsets_d), BATCH_SIZE*sizeof(int)));
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].outs_d), BATCH_SIZE*sizeof(int)*size_outs));
#else
    cudaSafeCall(cudaMalloc((void**)&((*s)[i].outs_d), BATCH_SIZE*sizeof(int)*2*tile_size));
#endif
    size_t free, total;
    cudaMemGetInfo(&free,&total);
    printf("%d MB free of total %d MB\n",free/1024/1024,total/1024/1024);
  }

  // set print buffer size (debug only)
  printf("NOTE: print buffer size is made larger\n");
  cudaSafeCall(cudaDeviceSetLimit(cudaLimitPrintfFifoSize, (1<<20)*50));
}

void GPU_close(std::vector<GPU_storage> *s, int num_threads){
  for(int i = 0; i < num_threads; ++i){
    cudaSafeCall(cudaFree((void*)((*s)[i].ref_seqs_d))); 
    cudaSafeCall(cudaFree((void*)((*s)[i].query_seqs_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].ref_lens_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].query_lens_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].ref_poss_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].query_poss_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].reverses_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].firsts_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].outs_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].matrices_d)));
#ifdef STREAM
    cudaSafeCall(cudaStreamDestroy((*s)[i].stream->stream));
    free((*s)[i].stream);
#endif
#ifdef GASAL
    cudaSafeCall(cudaFree((void*)((*s)[i].packed_ref_seqs_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].packed_query_seqs_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].ref_offsets_d)));
    cudaSafeCall(cudaFree((void*)((*s)[i].query_offsets_d)));
#endif
  }
}

#endif // GPU
