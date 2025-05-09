CC = g++
CFLAGS = -std=c++11 -Wall

.PHONY: all clean plot cls

all: main_q2

main_q2: main_q2.o ../source.o
	$(CC) $(CFLAGS) -o main_q2 main_q2.o ../source.o

main_q2.o: main_q2.cpp ../header.h
	$(CC) $(CFLAGS) -c main_q2.cpp

../source.o: ../source.cpp ../header.h
	$(CC) $(CFLAGS) -c ../source.cpp -o ../source.o

plot:
	./main_q2
	python3 plot_q2.py

clean:
	rm -f *.o main_q2 ../source.o

cls:
	rm -f *.o main_q2 ../source.o output.csv plots/*.png
	clear

