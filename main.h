#ifndef SECOND_H
#define SECOND_H

#include <vector>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <initializer_list>
#include <numeric>
#include <fenv.h>
#include "bmpLibrary/EasyBMP.h"

using namespace std; 


const int IMGWIDTH = 1920;
const int IMGHEIGHT = 1080;

/*
 * Function: applyHomo
 * ----------------------------
 *   Applies the homogenous matrix to a vector
 * 
 *   currPos: The a pointer to a vector
 *   homoMat: A pointer to a homgenous matrix
 * 
 */
void applyHomo(vector<double> *currPos, vector< vector<double> > *homoMat){
    double resultPos[4];

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            resultPos[i] += ((*currPos)[j] * (*homoMat)[j][i]);
        }
    }
    (*currPos)[0] = resultPos[0];
    (*currPos)[1] = resultPos[1];
    (*currPos)[2] = resultPos[2];
    (*currPos)[3] = resultPos[3];
}

struct position{
    double x;
    double y;
    double z;
    double w;
};

struct colour{
    double r;
    double g;
    double b;
};

class light{
private:
    struct position pos;
    struct position direction;
    struct colour col;
    double lightIntesnsity;
    string lightType;
public:
    vector<double> getPosition(){
        vector<double> position;
        position.push_back(pos.x);
        position.push_back(pos.y);
        position.push_back(pos.z);
        return(position);
    }

    vector<double> getDirection(){
        vector<double> direct;
        direct.push_back(direction.x);
        direct.push_back(direction.y);
        direct.push_back(direction.z);
        return(direct);
    }

    vector<double> getRGB(){
        vector<double> RGB;
        RGB.push_back(col.r);
        RGB.push_back(col.g);
        RGB.push_back(col.b);
        return(RGB);
    }
    void applyHomoLight( vector< vector<double> >  *homoMat){
        vector<double> currPos = getPosition();
        applyHomo(&currPos, homoMat);
        pos.x = currPos[0];
        pos.y = currPos[1];
        pos.z = currPos[2];
        pos.w = currPos[3];
    }
    void setPos(double posX, double posY, double posZ, double posW){
        pos.x = posX;
        pos.y = posY;
        pos.z = posZ;
        pos.w = posW;
    }
    double getLightIntensity(){
        return lightIntesnsity;
    }

    string getLightType(){
        return (lightType);
    }

    light(double posx, double posy, double posz, double directx, double directy, double directz, double r, double g, double b, double intensity, string tp);
};

light::light(double posx, double posy, double posz, double directx, double directy, double directz, double r, double g, double b, double intensity, string tp){
    pos.x = posx;
    pos.y = posy;
    pos.z = posz;
    direction.x = directx;
    direction.y = directy;
    direction.z = directz;
    col.r = r;
    col.g = g;
    col.b = b;
    lightIntesnsity = intensity;
    lightType = tp;
}


class vertex {
private:
    struct position pos;
    struct position winPos;
    struct position norm;
    struct colour col;
    struct colour specCol;
    struct colour diff;
    struct colour spec;
    int numbFaces = 0;
public:
    void setNormal(double xp, double yp, double zp){
        norm.x = xp;
        norm.y = yp;
        norm.z = zp;
    }
    
    void setNormal(vector<double> norms){
        norm.x = norms[0];
        norm.y = norms[1];
        norm.z = norms[2];
    }
    
    void incrementFaceCount(){
        numbFaces++;
    }

    void addDiffuse(double r, double g, double b){
        diff.r += r;
        diff.g += g;
        diff.b += b;
    }

    void setSpec(double r, double g, double b){
        spec.r = r;
        spec.g = g;
        spec.b = b;
    }

    void averageNormals(){
        norm.x = norm.x / (double) numbFaces;
        norm.y = norm.y / (double) numbFaces;
        norm.z = norm.z / (double) numbFaces;

        double length= sqrt(pow(norm.x, 2.0) + pow(norm.y, 2.0) + pow(norm.z, 2.0));
        norm.x /= length;
        norm.y /= length;
        norm.z /= length;
    }

    vector<double> getColour(){
        vector<double> colour;
        colour.push_back(col.r);
        colour.push_back(col.g);
        colour.push_back(col.b);
        return(colour);
    }

    vector<double> getDiff(){
        vector<double> colour;
        colour.push_back(diff.r);
        colour.push_back(diff.g);
        colour.push_back(diff.b);
        return(colour);
    }

    vector<double> getSpecColour(){
        vector<double> colour;
        colour.push_back(specCol.r);
        colour.push_back(specCol.g);
        colour.push_back(specCol.b);
        return(colour);
    }

    vector<double> getSpecAndDiff(){
        vector<double> colour;
        double red, green, blue;
        if (spec.r + diff.r > 255.0){
            red = 255.0;
        } else{
            red = spec.r + diff.r;
        }
        if (spec.g + diff.g > 255.0){
            green = 255.0;
        } else{
            green = spec.g + diff.g;
        }
        if (spec.b + diff.b > 255.0){
            blue = 255.0;
        } else{
            blue = spec.b + diff.b;
        }
        colour.push_back(red);
        colour.push_back(green);
        colour.push_back(blue);
        return(colour);
    }

    vector<double> getSpec(){
        vector<double> colour;
        colour.push_back(spec.r);
        colour.push_back(spec.g);
        colour.push_back(spec.b);
        return(colour);
    }

    vector<double> getNormal(){
        vector<double> normal;
        normal.push_back(norm.x);
        normal.push_back(norm.y);
        normal.push_back(norm.z);
        return(normal);
    }
    
    vector<double> getPosMat(){
        vector<double> position;
        position.push_back(pos.x);
        position.push_back(pos.y);
        position.push_back(pos.z);
        position.push_back(pos.w);
        return(position);
    }

    vector<double> getWinPosMat(){
        vector<double> position;
        position.push_back(winPos.x);
        position.push_back(winPos.y);
        position.push_back(winPos.z);
        position.push_back(winPos.w);
        return(position);
    }
    vector<int> getWinPosMatInt(){
        vector<int> position;
        position.push_back((int)(winPos.x + 0.5));
        position.push_back((int)(winPos.y + 0.5));
        position.push_back((int)(winPos.z + 0.5));
        position.push_back((int)(winPos.w + 0.5));
        return(position);
    }

    void applyHomoVert( vector< vector<double> >  *homoMat){
        vector<double> currPos = getPosMat();
        applyHomo(&currPos, homoMat);
        pos.x = currPos[0];
        pos.y = currPos[1];
        pos.z = currPos[2];
        pos.w = currPos[3];
    }

    void applyHomoVertWindow( vector< vector<double> >  *homoMat){
        vector<double> currPos = getWinPosMat();
        if (currPos[0] == NULL && currPos[1] == NULL && currPos[2] == NULL){
            currPos = getPosMat();
        }
        applyHomo(&currPos, homoMat);
        winPos.x = currPos[0];
        winPos.y = currPos[1];
        winPos.z = currPos[2];
        winPos.w = currPos[3];
    }

    void setPos(double posX, double posY, double posZ, double posW){
        pos.x = posX;
        pos.y = posY;
        pos.z = posZ;
        pos.w = posW;
    }

    void setWinPos(double posX, double posY, double posZ, double posW){
        winPos.x = posX;
        winPos.y = posY;
        winPos.z = posZ;
        winPos.w = posW;
    }
    vertex(double x, double y, double z, double r, double g, double b, double specr, double specg, double specb);
};

vertex::vertex(double x, double y, double z, double r, double g, double b, double specr, double specg, double specb){
    pos.x = x;
    pos.y = y;
    pos.z = z;
    pos.w = 1.0;
    norm.x = NULL;
    norm.y = NULL;
    norm.z = NULL;
    if (r > 255.0){
        col.r = 255.0;
    }else{
        col.r = r;
    }
    if (g > 255.0){
        col.g = 255.0;
    }else{
        col.g = g;
    }
    if (b > 255.0){
        col.b = 255.0;
    }else{
        col.b = b;
    }

    if (specr > 255.0){
        specCol.r = 255.0;
    }else{
        specCol.r = specr;
    }
    if (specg > 255.0){
        specCol.g = 255.0;
    }else{
        specCol.g = specg;
    }
    if (specb > 255.0){
        specCol.b = 255.0;
    }else{
        specCol.b = specb;
    }
    diff.r = 0.0;
    diff.g = 0.0;
    diff.b = 0.0;
    spec.r = 0.0;
    spec.g = 0.0;
    spec.b = 0.0;
    winPos.x = NULL;
    winPos.y = NULL;
    winPos.z = NULL;
}

class face {
private:
    struct position norm;
    struct colour col;
public:
    int v1;
    int v2;
    int v3;
    void setNormal(double xp, double yp, double zp){
        norm.x = xp;
        norm.y = yp;
        norm.z = zp;
    }

    void setNormal(vector<double> norms){
        norm.x = norms[0];
        norm.y = norms[1];
        norm.z = norms[2];
    }

    vector<double> getNormal(){
        vector<double> normal;
        normal.push_back(norm.x);
        normal.push_back(norm.y);
        normal.push_back(norm.z);
        return(normal);
    }
    face(int v1p, int v2p, int v3p);
};

face::face(int v1p, int v2p, int v3p){
    v1 = v1p;
    v2 = v2p;
    v3 = v3p;
}

class camera{
private:
    struct position pos;
public:
    vector<double> getPosMat(){
        vector<double> position;
        position.push_back(pos.x);
        position.push_back(pos.y);
        position.push_back(pos.z);
        position.push_back(pos.w);
        return(position);
    }
    camera(double posx, double posy, double posz);
};
camera::camera(double posx, double posy, double posz){
    pos.x = posx;
    pos.y = posy;
    pos.z = posz;
    pos.w = 1.0;
};

struct pixel{
    int r;
    int g;
    int b;
    float z;
    int faceNumb = -1;
};

/*
 * Function: dotProduct
 * ----------------------------
 *   Calculates the dot product
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: The dot product of two vectors
 * 
 */
double dotProduct(vector<double> mat1, vector<double> mat2){
    double result;
    result = (mat1[0] * mat2[0]) + (mat1[1] * mat2[1]) + (mat1[2] * mat2[2]);
    return result;
}

/*
 * Function: getDistanceLength
 * ----------------------------
 *   Calculates the distance length between two vectors
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: The distance length between two vectors
 * 
 */
double getDistanceLength(vector<double> mat1, vector<double> mat2){
    double distance = sqrt(pow(fabs(mat1[0] - mat2[0]), 2.0) +
                pow(fabs(mat1[1] - mat2[1]), 2.0) +
                pow(fabs(mat1[2] - mat2[2]), 2.0));

    return distance;
}

/*
 * Function: getUnitLength
 * ----------------------------
 *   Calculates the unit length between two vectors
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: A vector of the unit length between two vectors
 * 
 */
vector<double> getUnitLength(vector<double> mat1, vector<double> mat2){
    double distance = getDistanceLength(mat1, mat2);

    vector<double> result;
    result.push_back((mat2[0] - mat1[0]) / distance);
    result.push_back((mat2[1] - mat1[1]) / distance);
    result.push_back((mat2[2] - mat1[2]) / distance);

    return result;
}

/*
 * Function: getDirection
 * ----------------------------
 *   Calculates the direction vector between two vectors
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: The direction vector between two vectors
 * 
 */
vector<double> getDirection(vector<double> mat1, vector<double> mat2){
    vector<double> result;
    result.push_back((mat2[0] - mat1[0]) );
    result.push_back((mat2[1] - mat1[1]) );
    result.push_back((mat2[2] - mat1[2]));

    return result;
}

/*
 * Function: calculateHalfwayVector
 * ----------------------------
 *   Calculates the halfway vector between two vectors
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: The halfway vector between two vectors
 * 
 */
vector<double> calculateHalfwayVector(vector<double> mat1, vector<double> mat2){
    vector<double> result;
    double dx = mat1[0] - mat2[0];
    double dy = mat1[1] - mat2[1];
    double dz = mat1[2] - mat2[2];
    double denominator = sqrt(pow(dx, 2.0) + pow(dy, 2.0) + pow(dz, 2.0));
    result.push_back((mat1[0] + mat2[0]) / denominator);
    result.push_back((mat1[1] + mat2[1]) / denominator);
    result.push_back((mat1[2] + mat2[2]) / denominator);
    double distance = sqrt(pow(result[0], 2.0) + pow(result[1], 2.0) + pow(result[2], 2.0));
    result[0] /= distance;
    result[1] /= distance;
    result[2] /= distance;

    return result;
}

/*
 * Function: crossProduct
 * ----------------------------
 *   Calculates the cross product between two vectors
 * 
 *   mat1: A vector
 *   mat2: A vector
 * 
 *   return: The cross product between two vectors
 * 
 */
vector<double> crossProduct(vector<double> mat1, vector<double> mat2){
    vector<double> result;
    result.push_back((mat1[1] * mat2[2]) - (mat1[2] * mat2[1]));
    result.push_back((mat1[2] * mat2[0]) - (mat1[0] * mat2[2]));
    result.push_back((mat1[0] * mat2[1]) - (mat1[1] * mat2[0]));
    return result;
}

/*
 * Function: createHomoLine
 * ----------------------------
 *   Calculates a single line for a homogenous matrix
 * 
 *   v1: A double value for the matrix
 *   v2: A double value for the matrix
 *   v3: A double value for the matrix
 *   v4: A double value for the matrix
 * 
 *   return: A single line for a homogenous matrix
 * 
 */
vector<double> createHomoLine(double v1, double v2, double v3, double v4){
    vector<double> homoMatLine;
    homoMatLine.push_back(v1);
    homoMatLine.push_back(v2);
    homoMatLine.push_back(v3);
    homoMatLine.push_back(v4);
    return homoMatLine;
}

/*
 * Function: createHomoVect
 * ----------------------------
 *   Creates a 4x4 homogenous matrix
 * 
 *   homoMat: A pointer to a 2D vector which will become the homogenous matrix
 *   v1: A double value for the matrix
 *   v2: A double value for the matrix
 *   v3: A double value for the matrix
 *   v4: A double value for the matrix
 *   v5: A double value for the matrix
 *   v6: A double value for the matrix
 *   v7: A double value for the matrix
 *   v8: A double value for the matrix
 *   v9: A double value for the matrix
 *   v10: A double value for the matrix
 *   v11: A double value for the matrix
 *   v12: A double value for the matrix
 *   v13: A double value for the matrix
 *   v14: A double value for the matrix
 *   v15: A double value for the matrix
 *   v16: A double value for the matrix
 * 
 */
void createHomoVect(vector< vector<double> > *homoMat,
        double v1, double v2, double v3, double v4,
        double v5, double v6, double v7, double v8,
        double v9, double v10, double v11, double v12,
        double v13, double v14, double v15, double v16){
    homoMat->push_back(createHomoLine(v1, v5, v9, v13));
    homoMat->push_back(createHomoLine(v2, v6, v10, v14));
    homoMat->push_back(createHomoLine(v3, v7, v11, v15));
    homoMat->push_back(createHomoLine(v4, v8, v12, v16));

}

/*
 * Function: multiplyHomo
 * ----------------------------
 *   Multiplies together two homogenous matrices
 * 
 *   homoMat1: A homogenous matrix
 *   homoMat2: A homogenous matrix
 * 
 *   return: A homogenous matrix
 * 
 */
vector< vector<double> > multiplyHomo(vector< vector<double> > homoMat1, vector< vector<double> > homoMat2){
    vector< vector<double> > result;
    vector<double> homoMatLine;
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            homoMatLine.push_back((homoMat1[0][j] * homoMat2[i][0]) + 
                    (homoMat1[1][j] * homoMat2[i][1]) +
                    (homoMat1[2][j] * homoMat2[i][2]) +
                    (homoMat1[3][j] * homoMat2[i][3]));
                }
        result.push_back(homoMatLine);        
        homoMatLine.clear();
    }
    return result;
}

/*
 * Function: createRotateMatrix
 * ----------------------------
 *   Calculates a rotation homogenous matrix
 * 
 *   size: The size of the angle to rotate in degrees
 *   axis: The axis that is going to be rotated
 * 
 *   return: A homogenous matrix for rotation
 * 
 */
vector< vector<double> > createRotateMatrix(double size, char axis){
    vector< vector<double> > result;
    
    switch(axis){
        case 'x':
            createHomoVect(&result, 
        1, 0, 0, 0,
        0, cos(size * M_PI / 180), -sin(size* M_PI / 180), 0,
        0, sin(size* M_PI / 180), cos(size* M_PI / 180), 0,
        0, 0, 0, 1);
            break;
        case 'y':
            createHomoVect(&result, 
        cos(size* M_PI / 180), 0, sin(size* M_PI / 180), 0,
        0, 1, 0, 0,
        -sin(size* M_PI / 180), 0, cos(size* M_PI / 180), 0,
        0, 0, 0, 1);
            break;
        case 'z':
            createHomoVect(&result, 
        cos(size* M_PI / 180), -sin(size* M_PI / 180), 0, 0,
        sin(size* M_PI / 180), cos(size* M_PI / 180), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1);
            break;
    }

    return result;
}

/*
 * Function: createTranslateMatrix
 * ----------------------------
 *   Calculates the translation homogenous matrix
 * 
 *   x: Amount to translate in X direction
 *   y: Amount to translate in Y direction
 *   z: Amount to translate in Z direction
 * 
 *   return: A translation homogenous matrix
 * 
 */
vector< vector<double> > createTranslateMatrix(double x, double y, double z){
    vector< vector<double> > result;

    createHomoVect(&result, 
        1, 0, 0, x,
        0, 1, 0, y,
        0, 0, 1, z,
        0, 0, 0, 1);
    
    return result;
}

/*
 * Function: createScaleMatrix
 * ----------------------------
 *   Calculates the scale homogenous matrix
 * 
 *   x: Amount to scale in X 
 *   y: Amount to scale in Y 
 *   z: Amount to scale in Z 
 * 
 *   return: A scale homogenous matrix
 * 
 */
vector< vector<double> > createScaleMatrix(double x, double y, double z){
    vector< vector<double> > result;

    createHomoVect(&result, 
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1);
    
    return result;
}

/*
 * Function: edgeFunction
 * ----------------------------
 *   Calculates the edge function for two given vertices and a given x and y point
 * 
 *   v1: A pointer to a vertex
 *   v2: A pointer to a vertex
 *   px: A points X location
 *   py: A points Y location
 * 
 *   return: The edge function
 * 
 */
double edgeFunction(vertex* v1, vertex* v2, double px, double py){
    vector<double> v1Pos = v1->getWinPosMat();
    vector<double> v2Pos = v2->getWinPosMat();
    return((px - v1Pos[0]) * (v2Pos[1] - v1Pos[1]) - (py - v1Pos[1]) * (v2Pos[0] - v1Pos[0]));
}

/*
 * Function: edgeFunction
 * ----------------------------
 *   Calculates the edge function for two given point vectors and a given x and y point
 * 
 *   v1Pos: A specific vector of a point
 *   v2Pos: A specific vector of a point
 *   px: A points X location
 *   py: A points Y location
 * 
 *   return: The calculated red, green and blue light intensity
 * 
 */
double edgeFunction(vector<int> v1Pos, vector<int> v2Pos, double px, double py){
    return(((double)px - (double)v1Pos[0]) * ((double)v2Pos[1] - (double)v1Pos[1]) - ((double)py - (double)v1Pos[1]) * ((double)v2Pos[0] - (double)v1Pos[0]));
}

/*
 * Function: edgeFunction
 * ----------------------------
 *   Calculates the edge function for three given vertices
 * 
 *   v1: A pointer to a vertex
 *   v2: A pointer to a vertex
 *   v3: A pointer to a vertex
 * 
 *   return: The edge function
 * 
 */
double edgeFunction(vertex* v1, vertex* v2, vertex* v3){
    vector<double> v1Pos = v1->getWinPosMat();
    vector<double> v2Pos = v2->getWinPosMat();
    vector<double> v3Pos = v3->getWinPosMat();
    return((v3Pos[0] - v1Pos[0]) * (v2Pos[1] - v1Pos[1]) - (v3Pos[1] - v1Pos[1]) * (v2Pos[0] - v1Pos[0]));
}

/*
 * Function: checkEdgeFunction
 * ----------------------------
 *   Calculates checks if a point is within a calculated edge function of two vertices
 * 
 *   v1: A pointer to a vertex
 *   v2: A pointer to a vertex
 *   px: A points X location
 *   py: A points Y location
 * 
 *   return: A boolean saying if a point is within a calculated edge function of two vertices
 * 
 */
bool checkEdgeFunction(vertex* v1, vertex* v2, double px, double py){
    return(edgeFunction(v1, v2, px, py) >= 0);
}

/*
 * Function: orderRemaingVerts
 * ----------------------------
 *   Calculates the order for vertices positions, based on their y coordinates (from smallest to largest)
 * 
 *   v1Pos: A points vector position
 *   v2Pos: A points vector position
 *   vecPos: A pointer to a vector containing the ordered list of positions
 * 
 *   return: An integer defining the triangle type, 0 = mixed, 1 = top flat, 2 = bottom flat
 * 
 */
int orderRemaingVerts(vector<double> v1Pos, vector<double> v2Pos, vector< vector<double> >* vecPos){
    double threshold = 0.02;
    if ((v1Pos[1] - (*vecPos)[0][1] < threshold)){
        //Top flat
        vecPos->push_back(v1Pos);
        vecPos->push_back(v2Pos);
        return 1;
    } else if ((v2Pos[1] - (*vecPos)[0][1] < threshold)){
        //Top flat
        vecPos->push_back(v2Pos);
        vecPos->push_back(v1Pos);
        return 1;
    } else{
        if((fabs(v1Pos[1] - v2Pos[1]) < threshold)){
            //Bottom flat
            vecPos->push_back(v1Pos);
            vecPos->push_back(v2Pos);
            return 2;
        }else if (v1Pos[1] > v2Pos[1]){
            //Mixed
            vecPos->push_back(v2Pos);
            vecPos->push_back(v1Pos);
            return 0;
        } else if(v1Pos[1] < v2Pos[1]){
            //Mixed
            vecPos->push_back(v1Pos);
            vecPos->push_back(v2Pos);
            return 0;
        }
    }
}

/*
 * Function: calculateLineEquations
 * ----------------------------
 *   Calculates the gradient and y intercept for both the right and left line of the triangle
 * 
 *   triangleType: The triangle type
 *   v4x: A pointer to an additional vertex x location, used if triangle is mixed
 *   tempVert: A pointer to a vector containing the location of a tempary vertex, used if triangle is mixed
 * 
 */
void calculateLineEquations(int triangleType, int* v4x, vector<int>* tempVert,
     double* leftLineGradient, double* leftYIntercept, 
     double* rightLineGradient, double* rightYIntercept, 
     vector< vector<double> > vecPos){

    switch(triangleType){
            case 0:
                //Mixed triangle  
                (*v4x) = (int)(((double)vecPos[0][0] + ((((double)vecPos[1][1] - (double)vecPos[0][1]) / ((double)vecPos[2][1] - (double)vecPos[0][1])) * ((double)vecPos[2][0] - (double)vecPos[0][0]))) + 0.5);
                tempVert->push_back(*v4x);
                tempVert->push_back(vecPos[1][1]);
                if ((*tempVert)[0] > vecPos[1][0]){
                    //right side is longer
                    (*leftLineGradient) = (((double)vecPos[0][1] - (double)vecPos[1][1]) / ((double)vecPos[0][0] - (double)vecPos[1][0]));
                    (*leftYIntercept) = ((double)vecPos[0][1] - ((*leftLineGradient) * (double)vecPos[0][0]));
                    (*rightLineGradient) = ((double)(vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*rightYIntercept) = ((double)vecPos[0][1] - ((*rightLineGradient) * (double)vecPos[0][0]));
                }else{
                    //Left side longer
                    (*leftLineGradient) = (((double)vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*leftYIntercept) = ((double)vecPos[0][1] - ((*leftLineGradient) * (double)vecPos[0][0]));
                    (*rightLineGradient) = (((double)vecPos[0][1] - (double)vecPos[1][1]) / ((double)vecPos[0][0] - (double)vecPos[1][0]));
                    (*rightYIntercept) = ((double)vecPos[0][1] - ((*rightLineGradient) * (double)vecPos[0][0]));
                }
                break;
            case 1:
                //Top flat
                if (vecPos[1][0] > vecPos[0][0]){
                    (*leftLineGradient) = (((double)vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*leftYIntercept) = ((double)vecPos[0][1] - ((*leftLineGradient) * (double)vecPos[0][0]));
                    (*rightLineGradient) = (((double)vecPos[1][1] - (double)vecPos[2][1]) / ((double)vecPos[1][0] - (double)vecPos[2][0]));
                    (*rightYIntercept) = ((double)vecPos[1][1] - ((*rightLineGradient) * (double)vecPos[1][0]));
                }else{
                    (*leftLineGradient) = (((double)vecPos[1][1] - (double)vecPos[2][1]) / ((double)vecPos[1][0] - (double)vecPos[2][0]));
                    (*leftYIntercept) = ((double)vecPos[1][1] - ((*leftLineGradient) * (double)vecPos[1][0]));
                    (*rightLineGradient) = (((double)vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*rightYIntercept) = ((double)vecPos[0][1] - ((*rightLineGradient) * (double)vecPos[0][0]));
                }
                break;
            case 2:
                //Bottom Flat
                if (vecPos[1][0] < vecPos[2][0]){
                    (*leftLineGradient) = (((double)vecPos[0][1] - (double)vecPos[1][1]) / ((double)vecPos[0][0] - (double)vecPos[1][0]));
                    (*leftYIntercept) = ((double)vecPos[0][1] - ((*leftLineGradient) * (double)vecPos[0][0]));
                    (*rightLineGradient) = (((double)vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*rightYIntercept) = ((double)vecPos[0][1] - ((*rightLineGradient) * (double)vecPos[0][0]));
                }else{
                    (*leftLineGradient) = (((double)vecPos[0][1] - (double)vecPos[2][1]) / ((double)vecPos[0][0] - (double)vecPos[2][0]));
                    (*leftYIntercept) = ((double)vecPos[0][1] - ((*leftLineGradient) * (double)vecPos[0][0]));
                    (*rightLineGradient) = (((double)vecPos[0][1] - (double)vecPos[1][1]) / ((double)vecPos[0][0] - (double)vecPos[1][0]));
                    (*rightYIntercept) = ((double)vecPos[0][1] - ((*rightLineGradient) * (double)vecPos[0][0]));
                }
                break;
            }
        }
#endif