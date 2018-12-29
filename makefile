cpp = Function.cpp Function.1.cpp
cppo = Function.o Function.1.o
libPath = lib
lib = lib
liba = $(lib).a
makeName = main
makeFile = main.cpp

main:
	cd $(libPath) && g++ -c $(cpp) && ar -r lib$(liba) $(cppo)
	rm $(libPath)/*.o
	g++ -o $(makeName) $(makeFile) -L $(libPath) -l $(lib)

clean:
	rm $(libPath)/*.a
	rm main

