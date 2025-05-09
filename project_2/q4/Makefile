CC = g++
CFLAGS = -std=c++11 -Wall

.PHONY: all clean plot cls

all: main_q4

main_q4: main_q4.o ../source.o
	$(CC) $(CFLAGS) -o main_q4 main_q4.o ../source.o

main_q4.o: main_q4.cpp ../header.h
	$(CC) $(CFLAGS) -c main_q4.cpp

../source.o: ../source.cpp ../header.h
	$(CC) $(CFLAGS) -c ../source.cpp -o ../source.o

plot:
	python3 analysis.py

clean:
	rm -f *.o main_q4 ../source.o 

cls:
	rm -f *.o main_q4 ../source.o outputs/*.csv outputs/*.txt plots/*.png
	clear

