default:
	mpif90 -shared src/meshReader.f90 -o lib/libmeshReader.so -fPIC -Jlib/

clean:
	rm -rf lib/*.a lib/*.so lib/*.mod src/*.o
