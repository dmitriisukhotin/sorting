all: sort

sort: driver.o 
	g++ -o sort driver.o 

driver.o: driver.cpp 
	g++ -c driver.cpp

clean:
	rm sort *.o