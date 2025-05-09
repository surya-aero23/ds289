CC = g++
CFLAGS = -std=c++11 -Wall

.PHONY: all clean plot cls

all: main_q3

main_q3: main_q3.o ../source.o
	$(CC) $(CFLAGS) -o main_q3 main_q3.o ../source.o

main_q3.o: main_q3.cpp ../header.h
	$(CC) $(CFLAGS) -c main_q3.cpp

../source.o: ../source.cpp ../header.h
	$(CC) $(CFLAGS) -c ../source.cpp -o ../source.o

plot:
	./main_q3
	python3 plot_q3.py

clean:
	rm -f *.o main_q3 ../source.o

cls:
	rm -f *.o main_q3 ../source.o *.csv plots/*.png
	clear
