CPP= g++
CPPFLAGS = -O3 -std=c++11
SRCS_cpp = select.cpp
MAIN = selectFasta


all:    $(MAIN)
	@echo  RAPIFILT has been compiled
$(MAIN): $(OBJS1)
	$(CPP) $(CPPFLAGS) $(SRCS_cpp) -o $(MAIN) 

clean:
	$(RM) *.o *~ $(MAIN)

