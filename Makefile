default:
	mpif90 -c src/meshReader.f90 -o src/meshReader.o -Jsrc/
	ar crv lib/libmeshReader.a src/*.o

clean:
	rm -rf lib/*.a lib/*.so src/*.o src/*.mod
