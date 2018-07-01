
#CFLAGS = -O4 -g
CFLAGS = -O4

# linker flags, -pg enables use of gprof
#LFLAGS = -pg
LFLAGS =

#NCFLAGS = -O0 -G -Xptxas -v
NCFLAGS = -O3 -Xptxas -v

CFLAGS += $(options)

all: reference_guided 

ntcoding.o: ntcoding.h ntcoding.cpp
	g++ $(CFLAGS) -c ntcoding.cpp

fasta.o: fasta.h fasta.cpp 
	g++ $(CFLAGS) -c fasta.cpp

seed_pos_table.o: ntcoding.o seed_pos_table.h seed_pos_table.cpp
	g++ $(CFLAGS) -c seed_pos_table.cpp

reference_guided.o: reference_guided.cpp
	g++ -std=c++11 $(CFLAGS) -Wno-multichar -c reference_guided.cpp

chameleon.o: Chameleon.cpp
	g++ -std=c++11 $(CFLAGS) -Wno-multichar -o chameleon.o -c Chameleon.cpp

config_file.o: ConfigFile.cpp
	g++ -std=c++11 $(CFLAGS) -Wno-multichar -o config_file.o -c ConfigFile.cpp

reference_guided: ntcoding.o fasta.o seed_pos_table.o chameleon.o config_file.o reference_guided.o
	#g++ -std=c++11 $(CFLAGS) $(LFLAGS) -Wno-multichar -I/usr/local/include/ fasta.o ntcoding.o seed_pos_table.o chameleon.o config_file.o reference_guided.o cuda_host.cu gact.cpp align.cpp -o reference_guided -pthread 
	nvcc -std=c++11 $(NCFLAGS) $(LFLAGS) $(gpu_options) -arch=compute_35 -code=sm_35 -Xcompiler="-pthread -Wno-multichar" -I/usr/local/include/ fasta.o ntcoding.o seed_pos_table.o chameleon.o config_file.o reference_guided.cpp cuda_host.cu gact.cpp align.cpp -o reference_guided

clean:
	rm -rf *.o reference_guided 

