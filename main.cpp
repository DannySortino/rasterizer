#include "main.h"
#include "vertexLighting.cpp"
#include "fragmentLighting.cpp"


using namespace std; 
int type = 0; //0 = per vertex lighting, 1 = per fragment lighting - User parameter


/*
 * Function: readFile
 * ----------------------------
 *   Read the object file and store the face and vertex data
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: A pointer to a vector containing all faces
 *   file: The name and location of the object file to load
 * 
 *   return: an integer to state if there was any error in opening the file, 1 if there was an error, 0 if there was no error
 * 
 */
int readFile(vector<vertex>* verticesPtr, vector<face>* facesPtr, char* file){
    FILE* fp;
    char* characters = (char*)malloc(sizeof(char) * 200);

    if ((fp = fopen(file, "r")) == NULL){
        printf("ERROR! Cannot open file! \n");
        fclose(fp);
        return 1;
    }

    int v1, v2, v3;
    double x, y, z;
    
    while (fgets(characters, 199, fp) != NULL){
        switch(characters[0]){
            case 'v':
                if (sscanf(characters, "v %lf %lf %lf", &x, &y, &z) == 0){
                    printf("ERROR! Vertex line not defined correctly\n");
                    fclose(fp);
                    return 1;
                }
                verticesPtr->push_back(vertex(x, y, z, 100.0, 0.0, 0.0, 40.0, 40.0, 40.0));
                break;
            case 'f':
                if (sscanf(characters, "f %d %d %d",&v1, &v2, &v3) == 0){
                    printf("ERROR! Face line not defined correctly\n");
                    fclose(fp);
                    return 1;
                }
                facesPtr->push_back(face(v1 - 1, v2 - 1, v3 - 1));
                break;
        }
    }
    fclose(fp);
    return 0;
}

/*
 * Function: setBlackBackground
 * ----------------------------
 *   Sets the background of the image to black
 * 
 *   outputPtr: A pointer to the output BMP image
 * 
 */
void setBlackBackground(BMP *outputPtr){
    for (int i = 0; i < IMGWIDTH; i++){
        for (int j = 0; j < IMGHEIGHT; j++){
            (*outputPtr)(i, j)->Red = 0;
            (*outputPtr)(i, j)->Green = 0;
            (*outputPtr)(i, j)->Blue = 0;
        }
    }
}

/*
 * Function: applyHomoAllVerts
 * ----------------------------
 *   Applies the homogenous matrix to all vertices
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   homoMat: A pointer to a vector of size 4x4, corrosponding to a homogenous matrix
 * 
 */
void applyHomoAllVerts(vector<vertex>* verticesPtr, vector< vector<double> >* homoMat){
    for (int i = 0; i < verticesPtr->size(); i++){
        (*verticesPtr)[i].applyHomoVert(homoMat);
    }
}

/*
 * Function: applyHomoVertWindow
 * ----------------------------
 *   Applies the homogenous matrix to all vertices in the window space
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   homoMat: A pointer to a vector of size 4x4, corrosponding to a homogenous matrix
 * 
 */
void applyHomoVertWindow(vector<vertex>* verticesPtr, vector< vector<double> >* homoMat){
    for (int i = 0; i < verticesPtr->size(); i++){
        (*verticesPtr)[i].applyHomoVertWindow(homoMat);
    }
}

/*
 * Function: applyHomoAllLights
 * ----------------------------
 *   Applies the homogenous matrix to the light
 * 
 *   lightsPtr: A pointer to a vector of all lights
 *   homoMat: A pointer to a vector of size 4x4, corrosponding to a homogenous matrix
 * 
 */
void applyHomoAllLights(vector<light>* lightsPtr , vector< vector<double> >* homoMat){
    for (int i = 0; i < lightsPtr->size(); i++){
        (*lightsPtr)[i].applyHomoLight(homoMat);
    }
}

/*
 * Function: convertToWindowSpace
 * ----------------------------
 *   Converts the vertices from object space to window space
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 * 
 */
void convertToWindowSpace(vector<vertex>* verticesPtr, vector<light>* lightsPtr, camera cam){

    vector< vector<double> > homoMat, modelMat, viewMat, persMat;
    vector<double> camLoc = cam.getPosMat();
    //model matrix
    modelMat = createTranslateMatrix(0, 0, 2);

    //view matrix
    createHomoVect(&viewMat, 1, 0, 0, camLoc[0], 
                            0, 1, 0, camLoc[1], 
                            0 , 0, 1, camLoc[2],
                            0, 0, 0, 1);

    applyHomoVertWindow(verticesPtr, &viewMat);
    
    //perspective mat
    for (int i = 0; i < verticesPtr->size(); i++){
        vector<double> posMat = (*verticesPtr)[i].getWinPosMat();
        posMat[0] = posMat[0] / -posMat[2];
        posMat[1] = posMat[1] / -posMat[2];
        (*verticesPtr)[i].setWinPos((posMat[0] * 0.5 + 0.5) * IMGWIDTH, (posMat[1] * 0.5 + 0.5) * IMGHEIGHT, posMat[2], 1.0);
    }

    

}

/*
 * Function: createObject
 * ----------------------------
 *   Draws the framebuffer to the image
 * 
 *   outputPtr: A pointer to the output BMP image
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: A pointer to a vector containing all faces
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 * 
 */
void createObject(BMP *outputPtr, vector<vertex>* verticesPtr, vector<face>* facesPtr, vector<light>* lightsPtr, camera cam){
    vector< vector<pixel> > framebuffer(IMGWIDTH, vector<pixel>(IMGHEIGHT));
    convertToWindowSpace(verticesPtr, lightsPtr, cam);
    if (type == 0){
        rasterize(verticesPtr, facesPtr, &framebuffer);
    } else{
        fragRasterize(verticesPtr, facesPtr, &framebuffer, lightsPtr ,cam);
    }
    for (int x = 0; x < IMGWIDTH; x++){
        for (int y = 0; y < IMGHEIGHT; y++){
            (*outputPtr)(x, y)->Red = framebuffer[x][y].r;
            (*outputPtr)(x, y)->Green = framebuffer[x][y].g; 
            (*outputPtr)(x, y)->Blue = framebuffer[x][y].b; 
        }
    }
}

/*
 * Function: makeImage
 * ----------------------------
 *   Initalises the image file, calls other functions to fill the image, and then saves it
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: A pointer to a vector containing all faces
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 *   resultNum: The result number to generate the file
 * 
 */
void makeImage(vector<vertex>* verticesPtr, vector<face>* facesPtr, vector<light>* lightsPtr, camera cam, int resultNum){
    BMP output; 
    output.SetSize(IMGWIDTH, IMGHEIGHT);
    output.SetBitDepth(16);
    setBlackBackground(&output);
    createObject(&output, verticesPtr, facesPtr, lightsPtr, cam);
    char fileName[40];
    sprintf(fileName,"result/out_%d.bmp", resultNum);
    output.WriteToFile(fileName);
}

/*
 * Function: updateVetexNorm
 * ----------------------------
 *   Updates the vertex normal
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   vert: The vertex that is currently being updated
 *   normal: The normal that is being added to the vertex
 * 
 */
void updateVetexNorm(vertex* verticesPtr, vector<double> vert, vector<double> normal){
    verticesPtr->setNormal(vert[0] + normal[0], vert[1] + normal[1], vert[2] + normal[2]);
    verticesPtr->incrementFaceCount();    
}

/*
 * Function: calculateNormals
 * ----------------------------
 *   Calculates the normals for each vertex, based upon the faces it is connected to
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: A pointer to a vector containing all faces
 * 
 */
void calculateNormals(vector<vertex>* verticesPtr, vector<face>* facesPtr){
    vector<double> u, v, normal;
    for (int j = 0; j < facesPtr->size(); j++){
        vector<double> v1 = ((*verticesPtr)[(*facesPtr)[j].v1]).getPosMat();
        vector<double> v2 = ((*verticesPtr)[(*facesPtr)[j].v2]).getPosMat();
        vector<double> v3 = ((*verticesPtr)[(*facesPtr)[j].v3]).getPosMat();
        for(int i = 0; i < 4; i++){
            u.push_back(v2[i] - v1[i]);
            v.push_back(v3[i] - v1[i]);
        }
        normal = crossProduct(u, v);
        (*facesPtr)[j].setNormal(normal);

        vector<double> vert = ((*verticesPtr)[(*facesPtr)[j].v1]).getNormal();
        updateVetexNorm(&(*verticesPtr)[(*facesPtr)[j].v1], vert, normal);
        vert = ((*verticesPtr)[(*facesPtr)[j].v2]).getNormal();
        updateVetexNorm(&(*verticesPtr)[(*facesPtr)[j].v2], vert, normal);
        vert = ((*verticesPtr)[(*facesPtr)[j].v3]).getNormal();
        updateVetexNorm(&(*verticesPtr)[(*facesPtr)[j].v3], vert, normal);

        u.clear(); v.clear(); normal.clear();
    }
    for (int i = 0; i < verticesPtr->size(); i++){
        (*verticesPtr)[i].averageNormals();
    }
    
}

/*
 * Function: createLight
 * ----------------------------
 *   Creates a light object and adds it to the vector of lights
 * 
 *   lightsPtr: A pointer to a vector of all lights
 *   posx: The X position of the camera
 *   posy: The Y position of the camera
 *   posz: The Z position of the camera
 *   directx: The direction of the camera in the X direction
 *   directy: The direction of the camera in the Y direction
 *   directz: The direction of the camera in the Z direction
 *   r: The red colour of the light
 *   g: The green colour of the light
 *   b: The blue colour of the light
 *   lightInt: The Y position of the data element that is being read from
 *   tp: The string type of the light
 * 
 */
void createLight(vector<light>* lightsPtr, double posx, double posy, double posz,
     double directx, double directy, double directz,
      double r, double g, double b, double lightInt, string tp){
    lightsPtr->push_back(light(posx, posy, posz, directx, directy, directz, r, g, b,lightInt, tp));
}

/*
 * Function: main
 * ----------------------------
 *   Runs the main application, defining the variables needed, and calls the other functions to run the application
 * 
 *   argc: The number of parameters passed to the application
 *   argv: An array of the parameters passed to the application
 * 
 *   return: The state of execution from the main application (1 for error, 0 for ok)
 * 
 */
int main(int argc, char* argv[] ){
    vector<vertex> vertices; 
    vector<face> faces;
    vector<light> lights;
    camera cam(0.0, 0.0, 0.0);
    int imgNumb = -1;
    if (argc > 1){
        imgNumb = atoi(argv[1]);
        type = atoi(argv[2]);
    }
    //Load in obj file:
    char* fileToLoad = "StanfordBunny.obj";
    if (readFile(&vertices, &faces, fileToLoad) != 0){
        return 1;
    }
    printf("there are %d vertices\n", vertices.size());
    printf("there are %d faces\n", faces.size());
    createLight(&lights, -0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 255.0, 255.0, 255.0, 0.02, "point");
    vector< vector<double> > homoMat;
    homoMat = createRotateMatrix(imgNumb, 'y');
    applyHomoAllVerts(&vertices, &homoMat);
    homoMat = createTranslateMatrix(0, -0.1, 0.20);
    applyHomoAllVerts(&vertices, &homoMat);
    homoMat = createScaleMatrix(8.0, 8.0, 8.0);
    applyHomoAllVerts(&vertices, &homoMat);
    calculateNormals(&vertices, &faces);
    if (type == 0){
        calculateLighting(&vertices, &lights, cam);
    }
    makeImage(&vertices, &faces, &lights, cam, imgNumb);
    
    return 0;
}
