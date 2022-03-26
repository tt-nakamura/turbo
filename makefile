NTL = -lntl -lgmp
OBJ = TurboCode.o ConvolCode.o LDPCCode.o AddNoise.o GF2Xlib.o perm.o indexx.o

example: example.o $(OBJ)
	g++ example.o $(OBJ) $(NTL)
fig1: fig1.o $(OBJ)
	g++ fig1.o $(OBJ) $(NTL)
table1: table1.o $(OBJ)
	g++ table1.o $(OBJ) $(NTL)