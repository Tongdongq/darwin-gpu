
#include <stdio.h>

#ifndef CUDA_HEADER
#define CUDA_HEADER



#ifndef NOCUDA
__constant__ int _BATCH_SIZE;
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

typedef int AlnOp;
enum AlnOperands {ZERO_OP, DELETE_OP, INSERT_OP, MATCH_OP};
enum states {Z, D, I, M};

#define INF (1 << 30)

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
    const char reverse = reverses_d[tid];
    const char first = firsts_d[tid];
    const int row_len = _tile_size + 1;
    int *out = outs_d + tid * 2 * _tile_size;

    if(ref_len == -1){
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
    int *dir_matrix = d_matrix_rd + (_tile_size + 1) * NTHREADS;
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
        i_matrix_rd[i*__X] = -INF;
        d_matrix_rd[i*__X] = -INF;
       
        h_matrix_wr[i*__X] = 0;
        m_matrix_wr[i*__X] = 0;
        i_matrix_wr[i*__X] = -INF;
        d_matrix_wr[i*__X] = -INF;
    }

    //initialize the operands int he direction matrix for the first row and first
    //column
    for (int i = 0; i < ref_len + 1; i++) {
        dir_matrix[(i*row_len)*__X] = ZERO_OP;
    }

    for (int j = 0; j < query_len + 1; j++) {
        dir_matrix[j*__X] = ZERO_OP;
    }


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

        //j - row number; i - column number
        for (int j = 1; j < query_len + 1; j++) {
            //int ref_nt = (reverse) ? NtChar2Int(ref_seq[ref_len-i]) : NtChar2Int(ref_seq[i-1]);
            //int query_nt = (reverse) ? NtChar2Int(query_seq[query_len-j]) : NtChar2Int(query_seq[j-1]);
            // reverse indicates the direction of the alignment
            // 1: towards position = 0, 0: towards position = length
            int ref_nt = (reverse) ? ref_seq[(i-1)*__Y] : ref_seq[(ref_len-i)*__Y];
            int query_nt = (reverse) ? query_seq[(j-1)*__Y] : query_seq[(query_len-j)*__Y];
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
            tmp += (ins_open >= ins_extend) ? (2 << INSERT_OP) : 0;
            tmp += (del_open >= del_extend) ? (2 << DELETE_OP) : 0;
            (dir_matrix)[(i*row_len+j)*__X] = tmp;//*/
            

            /*if (h_matrix_wr[j] >= max_score) {
                max_score = h_matrix_wr[j];
                max_i = i;
                max_j = j;
            }// 4% speedup*/
            if (tmp2 >= max_score) {
                max_score = tmp2;
                max_i = i;
                max_j = j;
            }//*/

            /*if ((i == ref_pos) && (j == query_pos)) {
                pos_score = h_matrix_wr[j*__X];
            }//*/

        }
    }

    //printf("T%d tile done, max score: %d, max_i: %d, max_j: %d\n", tid, max_score, max_i, max_j);

    int *BT_states = out;
    int i = 1;
    int i_curr = ref_pos, j_curr = query_pos;
    int i_steps = 0, j_steps = 0;

    if(first){
        i_curr = max_i;
        j_curr = max_j;
        BT_states[i++] = max_score;
        BT_states[i++] = i_curr;
        BT_states[i++] = j_curr;
    }else{
        //BT_states[i++] = pos_score;
        BT_states[i++] = h_matrix_wr[query_len*__X];
    }

    int state = dir_matrix[(i_curr*row_len+j_curr)*__X] % 4;

    while (state != Z) {
        if ((i_steps >= _early_terminate) || (j_steps >= _early_terminate)) { // || (i_steps - j_steps > 30) || (i_steps - j_steps < -30)) {
            break;
        }
        BT_states[i++] = state;
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
    BT_states[0] = i - 1;
	//printf("tb done, i_curr: %d, j_curr: %d, i_steps: %d, j_steps: %d\n", \
    i_curr, j_curr, i_steps, j_steps);

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
