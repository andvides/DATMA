CC = gcc
CPP= g++
CPPFLAGS = -O3
CFLAGS= -D__USE_BSD -g
LFLAGS = -lz -lboost_iostreams
SRCS_cpp = process_fastqFile.cpp process_gzFile.cpp process_illumina_to_fastq.cpp \
		process_sff_to_fastq.cpp generic.cpp rapifilt.cpp 
SRCS_c =sff.c
OBJS1 = $(SRCS_cpp:.cpp=.o) 
OBJS2 = $(SRCS_c:.c=.o)

MAIN = rapifilt

.PHONY: depend clean

all:    $(MAIN)
	@echo  RAPIFILT has been compiled
$(MAIN): $(OBJS1) $(OBJS2) 
	$(CPP) $(CPPFLAGS) -o $(MAIN) $(OBJS1) $(OBJS2) $(LFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES)  -c $<  -o $@

.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES)  -c $<  -o $@


depend: $(SRCS_cpp) $(SRCS_c)
	makedepend $(INCLUDES) $^
	
clean:
	$(RM) *.o *~ $(MAIN)