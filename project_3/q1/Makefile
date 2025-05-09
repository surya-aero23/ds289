CC = g++
CFLAGS = -std=c++11 -Wall

.PHONY: all clean plot cls

all: main_q1

main_q1: main_q1.o ../source.o
	$(CC) $(CFLAGS) -o main_q1 main_q1.o ../source.o

main_q1.o: main_q1.cpp ../header.h
	$(CC) $(CFLAGS) -c main_q1.cpp

../source.o: ../source.cpp ../header.h
	$(CC) $(CFLAGS) -c ../source.cpp -o ../source.o

run:
	./main_q1
	
plot:
	python3 analysis.py

clean:
	rm -f *.o main_q1 ../source.o

cls:
	rm -f *.o main_q1 ../source.o *.csv
	clear
