CC = g++
CFLAGS = -g
LIBS = -lm

all: clean raytracer

raytracer:
	$(CC) $(CFLAGS) -c -o raytracer.o raytracer.cpp
	$(CC) $(CFLAGS) -c -o util.o util.cpp
	$(CC) $(CFLAGS) -c -o light_source.o light_source.cpp
	$(CC) $(CFLAGS) -c -o scene_object.o scene_object.cpp
	$(CC) $(CFLAGS) -c -o bmp_io.o bmp_io.cpp
	$(CC) $(CFLAGS) -o raytracer \
	raytracer.o util.o light_source.o scene_object.o bmp_io.o $(LIBS)

clean:
	-rm -f core *.o *.bmp
	-rm raytracer
	



