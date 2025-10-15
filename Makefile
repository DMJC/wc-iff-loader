build:
	g++ wc_viewer.cpp -std=c++17 -O2 -lSDL2 -lGLEW -lGL `pkg-config --libs --cflags openvr` -o wc_viewer
clean:
	rm wc_viewer
