TARGET = dst2root

CC      = g++
COPT    = -O2
LD      = gcc
LDFLAGS = $(OPT)
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
CXXFLAGS   = $(COPT) -Wall -fPIC $(ROOTCFLAGS)
CXXLIBS    = $(ROOTLIBS)
AUXFLAGS   = -L/usr/lib64 -lstdc++
LIBS          =  -lXpm -lX11 ${ROOTLIBS} -lm
OBJS = dst2root.o

all:$(TARGET)

$(TARGET): $(OBJS)
	@echo "======= Linking $(TARGET)"
	$(LD) -o $(TARGET) -O $(OBJS) $(LIBS) $(LDFLAGS) $(AUXFLAGS)
	@echo "======= Cleaning ..."
	rm -f $(OBJS)
	@echo "======= Done! "

clean:
	@echo "======= Cleaning everything."
	rm -f $(TARGET) $(OBJS)

.C.o:
	@echo "======= Generating Object ..."
	$(CC) $(CXXFLAGS) -c $<

dst2root.o: dst2root.C constant.h
