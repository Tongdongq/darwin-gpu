
all: reference_guided 

ntcoding.o: ntcoding.h ntcoding.cpp
	g++ -O4  -c ntcoding.cpp -g

fasta.o: fasta.h fasta.cpp 
	g++ -O4  -c fasta.cpp -g

seed_pos_table.o: ntcoding.o seed_pos_table.h seed_pos_table.cpp
	g++ -O4  -c seed_pos_table.cpp -g

reference_guided.o: reference_guided.cpp
	g++ -std=c++11 -O4 -Wno-multichar -c reference_guided.cpp -g

chameleon.o: Chameleon.cpp
	g++ -std=c++11 -O4 -Wno-multichar -o chameleon.o -c Chameleon.cpp -g

config_file.o: ConfigFile.cpp
	g++ -std=c++11 -O4 -Wno-multichar -o config_file.o -c ConfigFile.cpp -g

reference_guided: ntcoding.o fasta.o seed_pos_table.o chameleon.o config_file.o reference_guided.o 
	g++ -std=c++11 -O4 -Wno-multichar -pg -I/usr/local/include/ fasta.o ntcoding.o seed_pos_table.o chameleon.o config_file.o reference_guided.o -o reference_guided -pthread 

clean:
	rm -rf *.o reference_guided 

