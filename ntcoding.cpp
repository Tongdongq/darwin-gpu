
#include "ntcoding.h"
#include <stdlib.h>

int shape_pos[32];
int shape_size;
bool is_ignore_lower = false;
 
uint32_t NtChar2Int (char nt) {
    if (!is_ignore_lower) {
        nt = toupper(nt);
    }
    switch(nt) {
        case 'A': return A_NT;
        case 'C': return C_NT;
        case 'G': return G_NT;
        case 'T': return T_NT;
        case 'N': return N_NT;
        default: return N_NT;
    }
}

void SetIgnoreLower() {
    is_ignore_lower = true;
}

uint32_t KmerToIndex(std::string kmer) {
    uint32_t index = 0; 
    char nt;
    for (int i = 0; i < kmer.length(); i++) {
        nt = (kmer.at(i));
        index = (index << 2) + NtChar2Int(nt);
    }
    return index;
}

void GenerateShapePos (std::string shape) {
    shape_size = 0;
    for (int i = 0; i < shape.length(); i++) {
        if (shape[i] == '1') {
            shape_pos[shape_size++] = i;
        }
    }
}

uint32_t GetKmerIndexAtPos (char* sequence, uint32_t pos) {
    uint32_t kmer = 0;
    for (int i = 0; i < shape_size; i++) {
            uint32_t nt = NtChar2Int(sequence[pos+shape_pos[i]]);
            if (nt != N_NT) {
                kmer = (kmer << 2) + nt;
            }
            else {
                kmer = (1 << 31);
                break;
            }
    }
    return kmer;
}

uint32_t GetKmerIndexAtPos (std::string sequence, uint32_t pos) {
    uint32_t kmer = 0;
    for (int i = 0; i < shape_size; i++) {
            uint32_t nt = NtChar2Int(sequence[pos+shape_pos[i]]);
            if (nt != N_NT) {
                kmer = (kmer << 2) + nt;
            }
            else {
                kmer = (1 << 31);
                break;
            }
    }
    return kmer;
}

