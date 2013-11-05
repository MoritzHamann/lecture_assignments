CC = gcc

all: gradient_map hough_transform

gradient_map: gradient_map.c
	${CC} -O3 -o $@ $< -lm

hough_transform: hough_transform.c
	${CC} -O3 -o $@ $< -lm

clean:
	rm -f gradient_map
	rm -f hough_transform
