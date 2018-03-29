# C++ Rasterizer

This project was used for a coursework assignment, where we needed to load in an OBJ file, 
do all the required lighting calculations, and then rasterize, to output an image.

An example result image is shown:

![alt text](https://raw.githubusercontent.com/DannySortino/rasterizer/master/out_0.bmp)

# To compile

To compile, simply run the following command:

gcc -o -Wall main main.c

# To Run

To run a per vetex lighting of the bunny, simply execute:

./run.bat

To create a per fragment lighting, execute as following (where {image number} is an integer
representing the name of the output file and rotation of the bunny in degrees):

./main {image number} 1
