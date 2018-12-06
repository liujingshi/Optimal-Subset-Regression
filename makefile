cpp = Hello.cpp
cppo = Hello.o
libPath = lib
lib = lib
liba = $(lib).a
makeName = main
makeFile = main.cpp

main:
	cd $(libPath) && g++ -c $(cpp) && ar -r lib$(liba) $(cppo)
	rm  $(libPath)/*.o
	g++ -o $(makeName) $(makeFile) -L $(libPath) -l $(lib)

