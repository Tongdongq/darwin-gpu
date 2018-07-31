
/*
Copyright 2018 Yatish Turakhia, Gill Bejerano and William Dally
Copyright 2018 Tong Dong Qiu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdio.h>

#ifndef CUDA_HEADER
#define CUDA_HEADER


struct CUDA_Stream_Holder {
    cudaStream_t stream;  
};


#ifndef NOCUDA
__constant__ int _tile_size;
__constant__ int _tile_overlap;
__constant__ int _gap_open;
__constant__ int _gap_extend;
__constant__ int _match;
__constant__ int _mismatch;
__constant__ int _early_terminate;

#define NTHREADS (gridDim.x*blockDim.x)

#ifdef COALESCE_MATRICES
    #define __X NTHREADS
#else
    #define __X 1
#endif

#ifdef COALESCE_BASES
    #define __Y NTHREADS
#else
    #define __Y 1
#endif

// for GASAL: DELETE_OP means a query_base is 'aligned' to a gap

typedef int AlnOp;
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};

#define INF (1 << 29)

#define MAX_SEQ_LEN 324

__global__ void gasal_pack_kernel( \
    uint32_t* unpacked_query_batch, uint32_t* unpacked_target_batch, \
    uint32_t *packed_query_batch, uint32_t* packed_target_batch, \
    int query_batch_tasks_per_thread, int target_batch_tasks_per_thread, \
    uint32_t total_query_batch_regs, uint32_t total_target_batch_regs) {

    int32_t i;
    const int32_t tid = (blockIdx.x * blockDim.x) + threadIdx.x;//thread ID
    uint32_t n_threads = gridDim.x * blockDim.x;
if(tid==0){
    //printf("packing\n");
    for(i = 0; i < 50; ++i){
        //printf("%d\n", ((char*)unpacked_query_batch)[i]);
    }
}
    for (i = 0; i < query_batch_tasks_per_thread &&  (((i*n_threads)<<1) + (tid<<1) < total_query_batch_regs); ++i) {
        uint32_t *query_addr = &(unpacked_query_batch[(i*n_threads)<<1]);
        uint32_t reg1 = query_addr[(tid << 1)]; //load 4 bases of the query sequence from global memory
        uint32_t reg2 = query_addr[(tid << 1) + 1]; //load  another 4 bases
        //printf("T%d query_addr: %p\n", tid, query_addr+2*tid);
        uint32_t packed_reg = 0;
        packed_reg |= (reg1 & 15) << 28;        // ---
        packed_reg |= ((reg1 >> 8) & 15) << 24; //    |
        packed_reg |= ((reg1 >> 16) & 15) << 20;//    |
        packed_reg |= ((reg1 >> 24) & 15) << 16;//    |
        packed_reg |= (reg2 & 15) << 12;        //     > pack sequence
        packed_reg |= ((reg2 >> 8) & 15) << 8;  //    |
        packed_reg |= ((reg2 >> 16) & 15) << 4; //    |
        packed_reg |= ((reg2 >> 24) & 15);      //----
        uint32_t *packed_query_addr = &(packed_query_batch[i*n_threads]);
        packed_query_addr[tid] = packed_reg; //write 8 bases of packed query sequence to global memory
        //printf("T%d packed_query_addr: %p\n", tid, packed_query_addr+tid);
    }
/*if(tid==1){
    for(i = 0; i < 80; ++i){
        printf("%08x ", packed_query_batch[i]);
    }printf("\n");
}//*/
    for (i = 0; i < target_batch_tasks_per_thread &&  (((i*n_threads)<<1) + (tid<<1)) < total_target_batch_regs; ++i) {
        uint32_t *target_addr = &(unpacked_target_batch[(i * n_threads)<<1]);
        uint32_t reg1 = target_addr[(tid << 1)]; //load 4 bases of the target sequence from global memory
        uint32_t reg2 = target_addr[(tid << 1) + 1]; //load  another 4 bases
        uint32_t packed_reg = 0;
        packed_reg |= (reg1 & 15) << 28;        // ---
        packed_reg |= ((reg1 >> 8) & 15) << 24; //    |
        packed_reg |= ((reg1 >> 16) & 15) << 20;//    |
        packed_reg |= ((reg1 >> 24) & 15) << 16;//    |
        packed_reg |= (reg2 & 15) << 12;        //     > pack sequence
        packed_reg |= ((reg2 >> 8) & 15) << 8;  //    |
        packed_reg |= ((reg2 >> 16) & 15) << 4; //    |
        packed_reg |= ((reg2 >> 24) & 15);      //----
        uint32_t *packed_target_addr = &(packed_target_batch[i * n_threads]);
        packed_target_addr[tid] = packed_reg; //write 8 bases of packed target sequence to global memory
    }

} // end gasal_pack_kernel()


#define FIND_MAX(curr, gidx) \
    maxXY_y = (maxHH < curr) ? gidx : maxXY_y;\
    maxHH = (maxHH < curr) ? curr : maxHH;

__global__ void gasal_local_kernel( \
    uint32_t *packed_query_batch, uint32_t *packed_target_batch, \
    const int32_t *query_batch_lens, const int32_t *target_batch_lens, \
    int32_t *query_offsets, int32_t *target_offsets, \
    const int *query_poss, const int *ref_poss,
    int *out,
    const char *firsts, char *dir_matrix
    /*, \
    int32_t *score, int32_t *query_batch_end, int32_t *target_batch_end*/) {
        int32_t i, j, k, m, l;
        int32_t e, z;
        int32_t maxHH = 0;//initialize the maximum score to zero
        int32_t prev_maxHH = 0;
        int32_t subScore;
        int32_t ridx, gidx;
        short3 HD;
        short3 initHD = make_short3(0, 0, 0);
        int32_t maxXY_x = 0;
        int32_t maxXY_y = 0;
        const uint32_t tid = (blockIdx.x * blockDim.x) + threadIdx.x;//thread ID
        int32_t query_len = query_batch_lens[tid];
        int32_t ref_len = target_batch_lens[tid];
        if(ref_len == -1){return;}
        const int query_pos = query_poss[tid];
        const int ref_pos = ref_poss[tid];
        const int row_len = _tile_size + 1;
        out += (_tile_size * 2 * tid);
        //uint32_t packed_target_batch_idx = target_batch_offsets[tid] >> 3; //starting index of the target_batch sequence
        //uint32_t packed_query_batch_idx = query_batch_offsets[tid] >> 3;//starting index of the query_batch sequence
        packed_query_batch += query_offsets[tid] >> 3;
        packed_target_batch += target_offsets[tid] >> 3;
        //printf("T%d query offset: %d, %p, ref offset: %d, %p, query_len: %d, ref_len: %d\n", tid, query_offsets[tid], packed_query_batch, target_offsets[tid], packed_target_batch, query_len, ref_len);
        //printf("T%d query_pos: %d, ref_pos: %d\n", tid, query_pos, ref_pos);
        /*if(tid==14){
            printf("query: ");
            for(i = 0; i < query_len/8; ++i){
                printf("%08x ", packed_query_batch[i], packed_query_batch+i);
            }printf("\nref: ");
            for(i = 0; i < ref_len/8; ++i){
                printf("%08x ", packed_target_batch[i], packed_target_batch+i);
            }printf("\n");
            printf("highest ref address: %p\n", packed_target_batch+ref_len/8);
        }//*/
        int pos_score;        // score of last calculated corner
        const char first = firsts[tid];
        uint32_t query_batch_regs = (query_len >> 3) + (query_len&7 ? 1 : 0);//number of 32-bit words holding query_batch sequence
        uint32_t target_batch_regs = (ref_len >> 3) + (ref_len&7 ? 1 : 0);//number of 32-bit words holding target_batch sequence
        //printf("T%d packed_query_batch: %p, query_regs: %d, target_regs: %d\n", tid, packed_query_batch, query_batch_regs, target_batch_regs);
//if(tid==0)printf("match: %d, mismatch: %d, open: %d, extend: %d\n", _match, _mismatch, _gap_open, _gap_extend);
        //-----arrays for saving intermediate values------
        short3 global[MAX_SEQ_LEN];
        int32_t h[9];
        int32_t f[9];
        int32_t p[9];
        int32_t d[9];
        //--------------------------------------------
        for (i = 0; i < MAX_SEQ_LEN; i++) {
            global[i] = initHD;
        }
        dir_matrix += tid;
        for (int i = 0; i < ref_len + 1; i++) {
            dir_matrix[(i*row_len)*__X] = ZERO_OP;
        }
        for (int j = 0; j < query_len + 1; j++) {
            dir_matrix[j*__X] = ZERO_OP;
        }//*/
        //dir_matrix += (_tile_size+3)*__X;
        for (i = 0; i < target_batch_regs; i++) { //target_batch sequence in rows
            for (m = 0; m < 9; m++) {
                    h[m] = 0;
                    f[m] = 0;
                    p[m] = 0;
                    d[m] = 0;
            }
            register uint32_t gpac = packed_target_batch[i];//load 8 packed bases from target_batch sequence
            gidx = i << 3;
            ridx = 0;
            for (j = 0; j < query_batch_regs; j++) { //query_batch sequence in columns
                register uint32_t rpac = packed_query_batch[j];//load 8 bases from query_batch sequence

                //--------------compute a tile of 8x8 cells-------------------
                    for (k = 28; k >= 0; k -= 4) {
                        uint32_t rbase = (rpac >> k) & 15;//get a base from query_batch sequence
//if(tid==0){printf("base: %d\n", rbase);}
                        //-----load intermediate values--------------
                        HD = global[ridx];
                        h[0] = HD.x;
                        e = HD.y;
                        z = HD.z;
                        //-------------------------------------------
                        int32_t prev_hm_diff = h[0] + _gap_open;
#pragma unroll 8
                        for (l = 28, m = 1; m < 9; l -= 4, m++) {
                            int ii = i*8+m-1;
                            int jj = j*8+7-k/4;
                            uint32_t gbase = (gpac >> l) & 15;//get a base from target_batch sequence
if(tid==2){
    //printf("gbase: %d, rbase: %d\n", gbase, rbase);
    if((ii == 5 || ii == 6) && (jj == 2 | jj == 1)){
        //printf("i: %d, j: %d\n", ii, jj);
    }
}
                            //DEV_GET_SUB_SCORE_LOCAL(subScore, rbase, gbase);//check equality of rbase and gbase
                            //int ins_open = h[m-1] + _gap_open;
                            int ins_open = z + _gap_open;
                            int ins_extend = e + _gap_extend;
                            //int del_open = h[m] + _gap_open;
                            int del_open = d[m] + _gap_open;
                            int del_extend = f[m] + _gap_extend;
                            AlnOp tmp = ZERO_OP;
                            int tmp2 = 0;
                            subScore = (rbase == gbase) ? _match : _mismatch;
                            int match = p[m] + subScore;
                            z = match;
                            d[m] = match;
                            if(match > tmp2){
                                tmp2 = match;
                                tmp = MATCH_OP;
                            }
                            e = max(ins_open, ins_extend);
                            if(e > tmp2){
                                tmp2 = e;
                                tmp = INSERT_OP;
                            }
                            f[m] = max(del_open, del_extend);
                            if(f[m] > tmp2){
                                tmp2 = f[m];
                                tmp = DELETE_OP;
                            }
#ifdef OLD
                            tmp += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
                            tmp += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
#endif
int hm = h[m];
                            h[m] = tmp2;

                            //f[m] = max(h[m] + _gap_open, f[m] + _gap_extend);//whether to introduce or extend a gap in query_batch sequence
                            //h[m] = max(h[m], f[m]);
                            //h[m] = max(h[m], 0);
                            //e = max(prev_hm_diff, e + _gap_extend);//whether to introduce or extend a gap in target_batch sequence
                            //e = max(h[m-1] + _gap_open, e + _gap_extend);
                            //h[m] = max(h[m], e);
                            //FIND_MAX(h[m], gidx + (m-1));//the current maximum score and corresponding end position on target_batch sequence
 
if(tid==1){
    //printf("T%d new max: %d, i: %d, j: %d\n", tid, h[m], ii, jj);
}
                            if(h[m] == maxHH){  
                                if(ii > maxXY_y){
                                    maxXY_y = ii;
                                    maxXY_x = jj;
                                    maxHH = h[m];
                                }
                                if(ii == maxXY_y){
                                    if(jj > maxXY_x){
                                        maxXY_y = ii;
                                        maxXY_x = jj;
                                        maxHH = h[m];
                                    }
                                }
                            }
                            if(h[m] > maxHH){
                                maxXY_y = ii;
                                maxXY_x = jj;
                                maxHH = h[m]; 
                            }
if(tid==0){
    //printf("score: %d, i: %d, j: %d, ref: %d, query: %d, f[m]: %d, e: %d, p[m]: %d\n", h[m], i*8+m-1, j*8+7-k/4, gbase, rbase, f[m], e, p[m]+subScore);
    //printf("maxHH: %d, f[m]: %d, e: %d, p[m]: %d\n", maxHH, f[m], e, p[m]+subScore);
    if(ii > 260 && jj > 80 && jj < 120){
        //printf("i: %d, j: %d, ref: %d, query: %d\n", ii, jj, gbase, rbase);
        //printf("X1 i: %d, j: %d, score: %d, ref: %d, query: %d, dir: %d\n", ii, jj, h[m], gbase, rbase, tmp);
        //printf("X i: %d, j: %d, score: %d, ref: %d, query: %d, f[m]: %d, e: %d, p[m]: %d, ins_open: %d, ins_extend: %d\n", \
            i*8+m-1, j*8+7-k/4, h[m], gbase, rbase, f[m], e, p[m]+subScore, ins_open, ins_extend);
    }
}
                            p[m] = h[m-1];
                            int idx = ii*row_len + jj;
                            idx *= __X;
                            dir_matrix[idx] = tmp;
if(tid==0){
    if(ii > 240 && ii < 270 && jj > 110 && jj < 150){
        //printf("XT%d i: %d, j: %d, score: %d, ref: %d, query: %d, dir: %d, m: %d, io: %d, ie: %d, do: %d, de: %d, dir: %p, h[m]: %d\n", \
        tid, ii, jj, h[m], gbase, rbase, tmp, match, ins_open, ins_extend, del_open, del_extend, dir_matrix+idx, hm);
    }
}
                            //dir_matrix[(ii*row_len+jj)*__X] = tmp;
                            if(ii == ref_pos-1 && jj == query_pos-1){
if(first == 0){
    //printf("T%d pos_score: %d, qpos: %d, rpos: %d\n", tid, h[m], query_pos, ref_pos);
}
                                pos_score = h[m];
                            }
                        }
                        //----------save intermediate values------------
                        HD.x = h[m-1];
                        HD.y = e;
                        HD.z = z;
                        global[ridx] = HD;
                        //---------------------------------------------
                        //maxXY_x = (prev_maxHH < maxHH) ? ridx : maxXY_x;//end position on query_batch sequence corresponding to current maximum score
                        //prev_maxHH = max(maxHH, prev_maxHH);
                        ridx++;
                    } // end of 8x8 tile
                //-------------------------------------------------------
            } // end for all query bases
        } // end for all target bases

        //score[tid] = maxHH;//copy the max score to the output array in the GPU mem
        //query_batch_end[tid] = maxXY_x;//copy the end position on query_batch sequence to the output array in the GPU mem
        //target_batch_end[tid] = maxXY_y;//copy the end position on target_batch sequence to the output array in the GPU mem

if(1){
        //printf("T%d tile done, max score: %d, max_i: %d, max_j: %d\n", tid, maxHH, maxXY_y, maxXY_x);
}

        i = 1;
        int i_curr = ref_pos-1, j_curr = query_pos-1;
        int i_steps = 0, j_steps = 0;
//if(tid==2)printf("X T%d curr i: %d, j: %d\n", tid, i_curr, j_curr);

        if(first){
            i_curr = maxXY_y;
            j_curr = maxXY_x;
            out[i++] = maxHH;
            out[i++] = i_curr+1;
            out[i++] = j_curr+1;
        }else{
            //out[i++] = 32767;
            out[i++] = pos_score;
            //printf("T%d non first score: %d, query_pos: %d, ref_pos: %d\n", tid, pos_score, query_pos, ref_pos);
        }
//printf("T%d dir_matrix: %p\n", tid, dir_matrix);
        int idx = i_curr*row_len + j_curr;
        idx *= __X;
        char state = dir_matrix[idx] % 4;
        //char state = dir_matrix[(i_curr*row_len+j_curr)*__X] % 4;

#ifndef OLD
        out[i++] = state;
if(tid==0)printf("state: %d, i_curr: %d, j_curr: %d, steps: %d %d, i: %d\n", state, i_curr+1, j_curr+1, i_steps, j_steps, i);
        while (state != Z) {
            if ((i_steps >= _early_terminate) || (j_steps >= _early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
                break;
            }
            if (state == M) {
                i_curr--;
                j_curr--;
                i_steps++;
                j_steps++;
            }
            else if (state == I) {
                i_curr--;
                i_steps++;
            }
            else if (state == D) {
                j_curr--;
                j_steps++;
            }
            int idx = (i_curr*row_len+j_curr)*__X;
            state = dir_matrix[idx] % 4;
            out[i++] = state;
if(tid==0)printf("state: %d, i_curr: %d, j_curr: %d, steps: %d %d, i: %d\n", state, i_curr+1, j_curr+1, i_steps, j_steps, i);
//if(maxHH==63&&maxXY_y==139&&maxXY_x==170)printf("X T%d state: %d, i: %d, j: %d, steps i: %d, j: %d, i: %d, %p\n", tid, state, i_curr, j_curr, i_steps, j_steps, i, dir_matrix+idx);
        };
        if(state == Z || i_steps >= _early_terminate || j_steps >= _early_terminate){
            i--;
        }//*/
#else
    while (state != Z) {
        if ((i_steps >= _early_terminate) || (j_steps >= _early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
            break;
        }
        out[i++] = state;
//if(tid==0)printf("state: %d, i_curr: %d, j_curr: %d, steps: %d %d, i: %d\n", state, i_curr, j_curr, i_steps, j_steps, i);
        if (state == M) {
            int idx = ((i_curr-1)*row_len+j_curr-1);
            idx *= __X;
            char t = dir_matrix[idx];
//if(tid==0)printf("state M: %d, i_curr: %d, j_curr: %d, %p\n", t, i_curr, j_curr, dir_matrix+idx);
            state = dir_matrix[idx] % 4;//*/
            //state = (dir_matrix[((i_curr-1)*row_len+j_curr-1)*__X] % 4);
            i_curr--;
            j_curr--;
            i_steps++;
            j_steps++;
        }
        else if (state == I) {
            int idx = (i_curr*row_len+j_curr)*__X;
            char t = dir_matrix[idx];
//if(tid==0)printf("state I: %d, i_curr: %d, j_curr: %d, %p\n", t, i_curr, j_curr, dir_matrix+idx);
            state = (dir_matrix[(i_curr*row_len+j_curr)*__X] & (2 << INSERT_OP)) ? M : I;
            i_curr--;
            i_steps++;
        }
        else if (state == D) {
            int idx = (i_curr*row_len+j_curr)*__X;
            char t = dir_matrix[idx];
//if(tid==0)printf("state D: %d, i_curr: %d, j_curr: %d, %p\n", t, i_curr, j_curr, dir_matrix+idx);
            state = (dir_matrix[(i_curr*row_len+j_curr)*__X] & (2 << DELETE_OP)) ? M : D;
            j_curr--;
            j_steps++;
        }
    };
#endif


        out[0] = i - 1;
    /*printf("T%d tb done, i_curr: %d, j_curr: %d, i_steps: %d, j_steps: %d\n", \
    tid, i_curr, j_curr, i_steps, j_steps);*/
    //printf("T%d has %d elements\n", tid, i-1);

        return;
} // end gasal_local_kernel()



__global__ void Align_Kernel(const char *ref_seqs_d, const char *query_seqs_d, \
    const int *ref_lens_d, const int *query_lens_d, \
    const int *ref_poss_d, const int *query_poss_d, \
    const char *reverses_d, const char *firsts_d, int *outs_d, int *matricess_d){

    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    //const int idx_b = blockIdx.x;
    //const int idx_t = threadIdx.x;
#ifdef COALESCE_BASES
    const char *ref_seq = ref_seqs_d + tid;
    const char *query_seq = query_seqs_d + tid;
#else
    const char *ref_seq = ref_seqs_d + tid * _tile_size;
    const char *query_seq = query_seqs_d + tid * _tile_size;
#endif
    const int ref_len = ref_lens_d[tid];
    const int query_len = query_lens_d[tid];
    const int ref_pos = ref_poss_d[tid];
    const int query_pos = query_poss_d[tid];
    //printf("T%d query_pos: %d, ref_pos: %d\n", tid, query_pos, ref_pos);
    //const char reverse = reverses_d[tid];
    const char first = firsts_d[tid];
    const int row_len = _tile_size + 1;
    int *out = outs_d + tid * 2 * _tile_size;

    if(ref_len == -1){
        out[0] = 0;
        return;
    }

    /*if(tid==2){
        for(int i = 0; i < ref_len; ++i){
            if(i % 10 == 0){printf(" ");}
            printf("%c", ref_seq[i]);
        }printf("\n");
        for(int i = 0; i < query_len; ++i){
            if(i % 10 == 0){printf(" ");}
            printf("%c", query_seq[i]);
        }printf("\n");
    }//*/
    //printf("T%d lens: %d %d, reverse: %d, first: %d\n", tid, ref_len, query_len, reverse, first);
#ifdef COALESCE_MATRICES
    int *h_matrix_wr = matricess_d + tid;
    int *m_matrix_wr = h_matrix_wr + (_tile_size + 1) * NTHREADS;
    int *i_matrix_wr = m_matrix_wr + (_tile_size + 1) * NTHREADS;
    int *d_matrix_wr = i_matrix_wr + (_tile_size + 1) * NTHREADS;
    int *h_matrix_rd = d_matrix_wr + (_tile_size + 1) * NTHREADS;
    int *m_matrix_rd = h_matrix_rd + (_tile_size + 1) * NTHREADS;
    int *i_matrix_rd = m_matrix_rd + (_tile_size + 1) * NTHREADS;
    int *d_matrix_rd = i_matrix_rd + (_tile_size + 1) * NTHREADS;
    char *dir_matrix = (char*)(matricess_d + (_tile_size + 1) * 8 * NTHREADS) + tid;
    //printf("T%d h_wr: %p, m_wr: %p, h_rd: %p\n", tid, h_matrix_wr, m_matrix_wr, h_matrix_rd);
#else
    int *h_matrix_wr = matricess_d + (_tile_size + 1) * (_tile_size + 9) * tid;
    int *m_matrix_wr = h_matrix_wr + (_tile_size + 1);
    int *i_matrix_wr = m_matrix_wr + (_tile_size + 1);
    int *d_matrix_wr = i_matrix_wr + (_tile_size + 1);
    int *h_matrix_rd = d_matrix_wr + (_tile_size + 1);
    int *m_matrix_rd = h_matrix_rd + (_tile_size + 1);
    int *i_matrix_rd = m_matrix_rd + (_tile_size + 1);
    int *d_matrix_rd = i_matrix_rd + (_tile_size + 1);
    int *dir_matrix = d_matrix_rd + (_tile_size + 1);
#endif
    for (int i = 0; i < query_len + 1; i++) {
        h_matrix_rd[i*__X] = 0;
        m_matrix_rd[i*__X] = 0;
        //i_matrix_rd[i*__X] = -INF;
        //d_matrix_rd[i*__X] = -INF;
        i_matrix_rd[i*__X] = 0;
        d_matrix_rd[i*__X] = 0;
       
        h_matrix_wr[i*__X] = 0;
        m_matrix_wr[i*__X] = 0;
        //i_matrix_wr[i*__X] = -INF;
        //d_matrix_wr[i*__X] = -INF;
        i_matrix_wr[i*__X] = 0;
        d_matrix_wr[i*__X] = 0;
    }

    //initialize the operands int he direction matrix for the first row and first
    //column
    for (int i = 0; i < ref_len + 1; i++) {
        dir_matrix[(i*row_len)*__X] = ZERO_OP;
    }

    for (int j = 0; j < query_len + 1; j++) {
        dir_matrix[j*__X] = ZERO_OP;
    }
/*if(tid==14){
    printf("query: ");
    for(int i = 0; i < query_len; ++i){
        if(i%8==0)printf(" ");
        printf("%d", query_seq[i*__Y]);
    }printf("\n");
    printf("ref: ");
    for(int i = 0; i < ref_len; ++i){
        if(i%8==0)printf(" ");
        printf("%d", ref_seq[i*__Y]);
    }printf("\n");
}*/
//printf("query_seq: %p\n", query_seq);
    int max_score = 0;
    int max_i = 0;
    int max_j = 0;

    for (int i = 1; i < ref_len + 1; i++) {
        for (int k = 1; k < _tile_size + 1; k++) {
            m_matrix_rd[k*__X] = m_matrix_wr[k*__X];
            h_matrix_rd[k*__X] = h_matrix_wr[k*__X];
            i_matrix_rd[k*__X] = i_matrix_wr[k*__X];
            d_matrix_rd[k*__X] = d_matrix_wr[k*__X];
        }
char ref_nt = ref_seq[(i-1)*__Y];
        //j - row number; i - column number
        for (int j = 1; j < query_len + 1; j++) {
            //int ref_nt = (reverse) ? NtChar2Int(ref_seq[ref_len-i]) : NtChar2Int(ref_seq[i-1]);
            //int query_nt = (reverse) ? NtChar2Int(query_seq[query_len-j]) : NtChar2Int(query_seq[j-1]);
            // reverse indicates the direction of the alignment
            // 1: towards position = 0, 0: towards position = length
            //char ref_nt = (reverse) ? ref_seq[(i-1)*__Y] : ref_seq[(ref_len-i)*__Y];
            //char query_nt = (reverse) ? query_seq[(j-1)*__Y] : query_seq[(query_len-j)*__Y];
            //char ref_nt = ref_seq[(i-1)*__Y];
            char query_nt = query_seq[(j-1)*__Y];
            int match = (query_nt == ref_nt) ? _match : _mismatch;

            //columnwise calculations
            // find out max value
            /*if (m_matrix_rd[j-1] > i_matrix_rd[j-1] && m_matrix_rd[j-1] > d_matrix_rd[j-1]) {
                m_matrix_wr[j] = m_matrix_rd[j-1] + match;
            } else if (i_matrix_rd[j-1] > d_matrix_rd[j-1]) {
                m_matrix_wr[j] = i_matrix_rd[j-1] + match;
            } else {
                m_matrix_wr[j] = d_matrix_rd[j-1] + match;
            }
            if (m_matrix_wr[j] < 0) {
                m_matrix_wr[j] = 0;
            }//*/
            int tmp2 = m_matrix_rd[(j-1)*__X];
            if(i_matrix_rd[(j-1)*__X] > tmp2){
                tmp2 = i_matrix_rd[(j-1)*__X];
            }
            if(d_matrix_rd[(j-1)*__X] > tmp2){
                tmp2 = d_matrix_rd[(j-1)*__X];
            }
            tmp2 += match;
            if(tmp2 < 0){
                tmp2 = 0;
            }
            m_matrix_wr[j*__X] = tmp2;

            int ins_open   = m_matrix_rd[j*__X] + _gap_open;
            int ins_extend = i_matrix_rd[j*__X] + _gap_extend;
            int del_open   = m_matrix_wr[(j-1)*__X] + _gap_open;
            int del_extend = d_matrix_wr[(j-1)*__X] + _gap_extend;

            i_matrix_wr[j*__X] = (ins_open > ins_extend) ? ins_open : ins_extend;

            d_matrix_wr[j*__X] = (del_open > del_extend) ? del_open : del_extend;

            /*int max1 = m_matrix_wr[j] > i_matrix_wr[j] ? m_matrix_wr[j] : i_matrix_wr[j];
            int max2 = d_matrix_wr[j] > 0 ? d_matrix_wr[j] : 0;
            h_matrix_wr[j] = max1 > max2 ? max1 : max2;//*/

            /*(dir_matrix)[i*row_len+j] = ((m_matrix_wr[j] >= i_matrix_wr[j]) ? \
              ((m_matrix_wr[j] >= d_matrix_wr[j]) ? MATCH_OP : DELETE_OP) : \
              ((i_matrix_wr[j] >= d_matrix_wr[j]) ? INSERT_OP : DELETE_OP));

            if ((m_matrix_wr[j] <= 0) && (i_matrix_wr[j] <= 0) && (d_matrix_wr[j] <= 0)) {
                (dir_matrix)[i*row_len+j] = ZERO_OP;
            }//*/

            // this code is slower on CPU, probably due to lower ILP
            // but it is 13-14% faster on GPU, probably due to less branching
            AlnOp tmp = ZERO_OP;
            tmp2 = 0;
            if(m_matrix_wr[j*__X] > tmp2){
              tmp2 = m_matrix_wr[j*__X];
              tmp = MATCH_OP;
            }
            if(i_matrix_wr[j*__X] > tmp2){
              tmp2 = i_matrix_wr[j*__X];
              tmp = INSERT_OP;
            }
            if(d_matrix_wr[j*__X] > tmp2){
              tmp2 = d_matrix_wr[j*__X];
              tmp = DELETE_OP;
            }
            /*(dir_matrix)[(i*row_len+j)*__X] = tmp;
            h_matrix_wr[j*__X] = tmp2;
            (dir_matrix)[(i*row_len+j)*__X] += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
            (dir_matrix)[(i*row_len+j)*__X] += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
            // 7% speedup*/
            h_matrix_wr[j*__X] = tmp2;


#ifdef STABLE
            tmp += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
            tmp += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
#endif
            (dir_matrix)[(i*row_len+j)*__X] = tmp;//*/
if(tid==0){
    if(i-1 > 260 && j-1 > 80 && j-1 < 120){
        //printf("XT%d i: %d, j: %d, score: %d, ref: %d, query: %d, dir: %d, io: %d, ie: %d, do: %d, de: %d\n", \
        tid, i-1, j-1, tmp2, ref_nt, query_nt, tmp, ins_open, ins_extend, del_open, del_extend);
    }
}
            /*if (h_matrix_wr[j] >= max_score) {
                max_score = h_matrix_wr[j];
                max_i = i;
                max_j = j;
            }// 4% speedup*/
if(tid==1){
    if(tmp2 >= max_score){
        //printf("T%d tmp2: %d, max_score: %d, i: %d, j: %d\n", tid, tmp2, max_score, i-1, j-1);
    }
}
            if (tmp2 >= max_score) {
                max_score = tmp2;
                max_i = i;
                max_j = j;
            }//*/
if(tid==23){
    if(i-1 < 1 && j-1 < 1){
        //printf("i: %d, j: %d, ref: %d, query: %d\n", i-1, j-1, ref_nt, query_nt);
        //printf("X i: %d, j: %d, score: %d, ref: %d, query: %d, m: %d, i: %d, d: %d\n", i-1, j-1, tmp2, ref_nt, query_nt, m_matrix_wr[j*__X], i_matrix_wr[j*__X], d_matrix_wr[j*__X]);
    }
}
if(tid==6){
    //if(i-1 > 300 && j-1 > 200){
        //printf("XT%d i: %d, j: %d, score: %d, ref: %d, query: %d, dir: %d, io: %d, ie: %d, do: %d, de: %d\n",\
     tid, i-1, j-1, tmp2, ref_nt, query_nt, tmp, ins_open, ins_extend, del_open, del_extend);
    //}
}
            if ((i == ref_pos) && (j == query_pos)) {
                //pos_score = h_matrix_wr[j*__X];
if(first == 0){
    //printf("T%d pos_score: %d, qpos: %d, rpos: %d\n", tid, h_matrix_wr[j*__X], query_pos, ref_pos);
}
            }//*/

        }
    }

    //printf("T%d tile done, max score: %d, max_i: %d, max_j: %d\n", tid, max_score, max_i-1, max_j-1);

    int *BT_states = out;
    int i = 1;
    int i_curr = ref_pos, j_curr = query_pos;
    int i_steps = 0, j_steps = 0;
//if(tid==2)printf("X T%d curr i: %d, j: %d\n", tid, i_curr-1, j_curr-1);

    if(first){
        i_curr = max_i;
        j_curr = max_j;
        BT_states[i++] = max_score;
        BT_states[i++] = i_curr;
        BT_states[i++] = j_curr;
    }else{
        //BT_states[i++] = pos_score;
        BT_states[i++] = h_matrix_wr[query_len*__X];
        //printf("T%d non first score: %d, query_pos: %d, ref_pos: %d\n", tid, h_matrix_wr[query_len*__X], query_pos, ref_pos);
    }

    char state = dir_matrix[(i_curr*row_len+j_curr)*__X] % 4;

#ifdef STABLE
    while (state != Z) {
        if ((i_steps >= _early_terminate) || (j_steps >= _early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
            break;
        }
        BT_states[i++] = state;
if(tid==0)printf("state: %d, i_curr: %d, j_curr: %d, steps: %d %d, i: %d\n", state, i_curr, j_curr, i_steps, j_steps, i);
        if (state == M) {
            state = (dir_matrix[((i_curr-1)*row_len+j_curr-1)*__X] % 4);
            i_curr--;
            j_curr--;
            i_steps++;
            j_steps++;
        }
        else if (state == I) {
            state = (dir_matrix[(i_curr*row_len+j_curr)*__X] & (2 << INSERT_OP)) ? M : I;
            i_curr--;
            i_steps++;
        }
        else if (state == D) {
            state = (dir_matrix[(i_curr*row_len+j_curr)*__X] & (2 << DELETE_OP)) ? M : D;
            j_curr--;
            j_steps++;
        }
    };
#else
    BT_states[i++] = state;
    while (state != Z) {
//if(tid==1&&query_pos==206)printf("X T%d state: %d, i: %d, j: %d, steps i: %d, j: %d, i: %d\n", tid, state, i_curr, j_curr, i_steps, j_steps, i);
        if ((i_steps >= _early_terminate) || (j_steps >= _early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
            break;
        }
        if (state == M) {
            i_curr--;
            j_curr--;
            i_steps++;
            j_steps++;
        }
        else if (state == I) {
            i_curr--;
            i_steps++;
        }
        else if (state == D) {
            j_curr--;
            j_steps++;
        }
        state = (dir_matrix[(i_curr*row_len+j_curr)*__X] % 4);
        BT_states[i++] = state;
    };
    if(state == Z || i_steps >= _early_terminate || j_steps >= _early_terminate){
        i--;
    }
#endif // STABLE
    BT_states[0] = i - 1;
    //printf("T%d tb done, i_curr: %d, j_curr: %d, i_steps: %d, j_steps: %d\n", \
    tid, i_curr-1, j_curr-1, i_steps, j_steps);

    /*std::queue<int> BT_states;
    int i_curr=ref_pos, j_curr=query_pos;
    int i_steps = 0, j_steps = 0;

    int open = 0;
    if (first) {
        i_curr = max_i;
        j_curr = max_j;
        BT_states.push(max_score);
        BT_states.push(i_curr);
        BT_states.push(j_curr);
    }
    else {
        BT_states.push(pos_score);
    }

    printf("tile done, first: %d, max_i: %d, max_j: %d, i_curr: %d, j_curr: %d, max_score: %d, pos_score: %d\n", \
    first, max_i, max_j, i_curr, j_curr, max_score, pos_score);

    int state = dir_matrix[i_curr][j_curr] % 4;

    while (state != Z) {
        if ((i_steps >= early_terminate) || (j_steps >= early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
            break;
        }
        BT_states.push(state);
        if (state == M) {
            state = (dir_matrix[i_curr-1][j_curr-1] % 4);
            i_curr--;
            j_curr--;
            i_steps++;
            j_steps++;
        }
        else if (state == I) {
            state = (dir_matrix[i_curr][j_curr] & (2 << INSERT_OP)) ? M : I;
            i_curr--;
            i_steps++;
        }
        else if (state == D) {
            state = (dir_matrix[i_curr][j_curr] & (2 << DELETE_OP)) ? M : D;
            j_curr--;
            j_steps++;
        }
    }; 

    printf("T%d tb done, i_curr: %d, j_curr: %d, i_steps: %d, j_steps: %d\n\n", \
    tid, i_curr, j_curr, i_steps, j_steps);//*/

    return;
}









// error checking code taken from: https://codeyarns.com/2011/03/02/how-to-do-error-checking-in-cuda/
#define cudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define cudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
  {
      if ( cudaSuccess != err )
      {
          //std::cerr << "cudaSafeCall() failed at " << file << ":" << line << ":" << cudaGetErrorString(err) << std::endl;
          printf("\ncudaSafeCall() failed at %s: %d: %s\n\n", file, line, cudaGetErrorString(err));
          exit( -1 );
      }
  }
inline void __cudaCheckError( const char *file, const int line )
  {
      cudaError err = cudaGetLastError();
      if ( cudaSuccess != err )
      {
          //std::cerr << "cudaCheckError() failed at " << file << ":" << line << ":" << cudaGetErrorString(err) << std::endl;
          printf("\ncudaSafeCall() failed at %s: %d: %s\n\n", file, line, cudaGetErrorString(err));
          exit( -1 );
      }
  }
// end of error checking code
#endif // NOCUDA

#endif
