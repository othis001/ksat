#lab3.make
# example make file for use in later labs

# first name a variable objects for all object files

objects = main.o Complex.o Smooshable.o


# name a variable sources for all source files

sources = main.cpp Complex.cpp Smooshable.cpp

# now give target as lab1 with objects as variable dependencies + command line

findEuler: $(objects) 
	g++ -o smoosher $(objects)

# how to build all object files from all dependent source files


$(objects): $(sources)
	g++ -c $(sources)

clean:
	rm $(objects) smoosher 

