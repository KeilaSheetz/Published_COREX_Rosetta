CC		= gcc
AR		= ar
ARFLAGS		= ruc
INCLUDES        = -I.
#CFLAGS		= -g -DDEBUG
CFLAGS		= $(INCLUDES) -O -c -O2
#CFLAGS		= -g
LFLAGS		=
LIBRARIES	= -lm -O
RM              = /bin/rm

EXE		= libasa.a EnsembleGeneratorMC
BINARY		= EnsembleGeneratorMC
OBJECTS		= EnsembleGeneratorMC.o libasa.a
ASA_OBJ		= libasa.o area_bugfix.o

all:	$(EXE)

$(EXE):	$(OBJS)

libasa.a: $(ASA_OBJ)
	$(AR) $(ARFLAGS) libasa.a $(ASA_OBJ)

EnsembleGeneratorMC: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o $(BINARY) $(LIBRARIES)

EnsembleGeneratorMC.o: EnsembleGeneratorMC.c asa_types.h Makefile
	$(CC) -c $(CFLAGS) EnsembleGeneratorMC.c

area_bugfix.o: area_bugfix.c libasa.h asa_types.h Makefile
	$(CC) -c $(CFLAGS) area_bugfix.c

libasa.o: libasa.c libasa.h asa_types.h radii.h Makefile
	$(CC) -c $(CFLAGS) libasa.c
