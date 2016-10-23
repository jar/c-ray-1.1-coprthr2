EXE = c-ray-mt

INCS = -I. -I/usr/local/browndeer/coprthr2/include
LIBS = -L/usr/local/browndeer/coprthr2/lib -lcoprthr -lcoprthrcc -lm

CC = gcc
CFLAGS = -O3
EDEFS = 

TARGET = $(EXE) device.e32

all: $(TARGET)

ifdef USESDL
INCS +=  -I/usr/local/include/SDL2
LIBS += -lSDL2 -lpthread
CFLAGS += -DUSESDL
endif

.PHONY: clean distclean
.SUFFIXES:
.SUFFIXES: .c .o .x .e32

$(EXE): main.o
	$(CC) -o $(EXE) main.o $(LIBS)

.c.e32:
	coprcc $(EDEFS) $< -o $@

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	rm -f *.o
	rm -f *.cbin.*

distclean: clean
	rm -f *.ppm
	rm -f $(TARGET) 
