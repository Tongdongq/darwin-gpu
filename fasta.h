#ifndef FASTA_H_
#define FASTA_H_

#include <string>
#include <vector>

#define SEQLINE_WRAP_LEN 70

void ParseFastaFile(std::string filename, 
                    std::vector<std::vector<std::string> >& descrips,
                    std::vector<std::string>& seqs,
                    std::vector<long long int>& lengths,
                    std::vector<long long int>& fileposs);

#endif // FASTA_H_
