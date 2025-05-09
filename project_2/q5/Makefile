CC = g++
CFLAGS = -std=c++11 -Wall

.PHONY: all clean plot cls

all: main_q5

main_q5: main_q5.o ../source.o
	$(CC) $(CFLAGS) -o main_q5 main_q5.o ../source.o

main_q5.o: main_q5.cpp ../header.h
	$(CC) $(CFLAGS) -c main_q5.cpp

../source.o: ../source.cpp ../header.h
	$(CC) $(CFLAGS) -c ../source.cpp -o ../source.o

plot:
	python3 analysis.py

run:
	./main_q5

clean:
	rm -f *.o main_q5 ../source.o

cls:
	rm -f *.o main_q5 ../source.o *.csv plots/*.png
	clear

