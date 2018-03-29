#include "main.h"

using namespace std; 


/*
 * Function: calculateLightIntensity
 * ----------------------------
 *   Calculates the intensity of a light for a given vertex
 * 
 *   vertexPt: A specific vertex
 *   lightPt: A specific light
 * 
 *   return: The calculated red, green and blue light intensity
 * 
 */
vector<double> calculateLightIntensity(vertex vertexPt, light lightPt){
    vector<double> vertexPosMat = vertexPt.getPosMat();
    vector<double> lightPosMat = lightPt.getPosition();
    double lightInt = lightPt.getLightIntensity();
    double distance = getDistanceLength(vertexPosMat, lightPosMat);
    double intensity = lightInt/pow(distance, 2.0);

    vector<double> lightCol = lightPt.getRGB();
    vector<double> RGB;

    RGB.push_back(lightCol[0] * intensity);
    RGB.push_back(lightCol[1] * intensity);
    RGB.push_back(lightCol[2] * intensity);

    return RGB;
}

/*
 * Function: addAmbient
 * ----------------------------
 *   Adds the ambient lighting to each vertex
 * 
 *   vertexPt: The number of parameters passed to the application
 *   r: The red colour of the ambient
 *   g: The green colour of the ambient
 *   b: The blue colour of the ambient
 * 
 */
void addAmbient(vertex* vertexPt, double r, double g, double b){
    vector<double> colour = vertexPt->getColour();
    vertexPt->addDiffuse(r * colour[0], g* colour[1], b* colour[2]);
    colour = vertexPt->getDiff();
}

/*
 * Function: calculateDiffuse
 * ----------------------------
 *   Calculates the diffuse for each vertex
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   lightsPtr: A pointer to a vector of all lights
 * 
 */
void calculateDiffuse(vector<vertex>* verticesPtr, vector<light>* lightsPtr){
    //for all verts
    for (int i = 0; i < verticesPtr->size(); i++){
        vector<double> lightIntensity (4, 0.0);
        for (int j = 0; j < lightsPtr->size(); j++){
            vector<double> tempLight = calculateLightIntensity((*verticesPtr)[i], (*lightsPtr)[j]);

            vector<double> vertNorm = (*verticesPtr)[i].getNormal();

            vector<double> vertexPosMat = (*verticesPtr)[i].getPosMat();
            vector<double> lightPosMat = (*lightsPtr)[j].getPosition();
            vector<double> lightUnitVec = getUnitLength(vertexPosMat, lightPosMat);

            double vertDot = dotProduct(vertNorm, lightUnitVec);

            if (vertDot < 0.0){
                vertDot = 0.0;
            }
            tempLight[0] = tempLight[0] * vertDot;
            tempLight[1] = tempLight[1] * vertDot;
            tempLight[2] = tempLight[2] * vertDot;

            //Used for vector summation
            transform(lightIntensity.begin(), lightIntensity.end(), tempLight.begin(), lightIntensity.begin(), plus<float>());
        }
        vector<double> colour = (*verticesPtr)[i].getColour();
        addAmbient(&((*verticesPtr)[i]), 0.3, 0.3, 0.3);
        (*verticesPtr)[i].addDiffuse((colour[0] * lightIntensity[0]), 
                                        (colour[1] * lightIntensity[1]), 
                                        (colour[2] * lightIntensity[2]));


        colour = (*verticesPtr)[i].getDiff();

    }
}

/*
 * Function: calculateSpecular
 * ----------------------------
 *   Calculates the specular for each vertex
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 *   m: The specular exponent
 * 
 */
void calculateSpecular(vector<vertex>* verticesPtr, vector<light>* lightsPtr, camera cam, double m){
    vector<double> cameraPos = cam.getPosMat();
    for (int i = 0; i < verticesPtr->size(); i++){
        vector<double> lightIntensity (4, 0.0);
        for (int j = 0; j < lightsPtr->size(); j++){

            vector<double> tempLight = calculateLightIntensity((*verticesPtr)[i], (*lightsPtr)[j]);
            vector<double> vertNorm = (*verticesPtr)[i].getNormal();
            vector<double> vertexPosMat = (*verticesPtr)[i].getPosMat();
            vector<double> lightPosMat = (*lightsPtr)[j].getPosition();
            vector<double> lightUnitVec = getUnitLength(vertexPosMat, lightPosMat);
            vector<double> camUnitVec = getUnitLength(vertexPosMat, cameraPos);

            double vertFaceLight = dotProduct(vertNorm, lightUnitVec);
                        
            //If normal dot (vector from point to light)
            if (vertFaceLight > 0.0){

                vector<double> halfwayVec = calculateHalfwayVector(lightUnitVec, camUnitVec);

                double vertDot = dotProduct(vertNorm, halfwayVec);

                if (vertDot < 0.0){
                    vertDot = 0.0;
                }
                double vertspec = pow(vertDot, m);
                tempLight[0] = tempLight[0] * vertspec;
                tempLight[1] = tempLight[1] * vertspec;
                tempLight[2] = tempLight[2] * vertspec;

            
                //Used for vector summation
                transform(lightIntensity.begin(), lightIntensity.end(), tempLight.begin(), lightIntensity.begin(), plus<float>());

            }
        vector<double> specColour = (*verticesPtr)[i].getSpecColour();
        (*verticesPtr)[i].setSpec((specColour[0] * lightIntensity[0]), 
                                        (specColour[1] * lightIntensity[1]), 
                                        (specColour[2] * lightIntensity[2]));
            
        }
    }
}

/*
 * Function: rasterize
 * ----------------------------
 *   Runs the main rasterising code, using a scanline approach.
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: An array of the parameters passed to the application
 *   framebuffer: A 2D vector containing the frame buffer data that will be output to the image
 * 
 */
void rasterize(vector<vertex>* verticesPtr, vector<face>* facesPtr, vector< vector<pixel> >* framebuffer){
    double near = 0.1;
    double far = 20.0;
    for (int i = 0; i < facesPtr->size(); i++){
        vertex* v1  = &(*verticesPtr)[(*facesPtr)[i].v1];
        vertex* v2  = &(*verticesPtr)[(*facesPtr)[i].v2];
        vertex* v3  = &(*verticesPtr)[(*facesPtr)[i].v3];

        vector<double> v1Pos = v1->getWinPosMat();
        vector<double> v2Pos = v2->getWinPosMat();
        vector<double> v3Pos = v3->getWinPosMat();

        if (v1Pos[2] < near || v2Pos[2] < near || v3Pos[2] < near){
            continue;
        }

        if (v1Pos[2] > far || v2Pos[2] > far || v3Pos[2] > far){
            continue;
        }

        double area = edgeFunction(v1, v2, v3);
        int triangleType;
        
        vector< vector<double> > vecPos;
        if (v1Pos[1] <= v2Pos[1] && v1Pos[1] <= v3Pos[1]){
            vecPos.push_back(v1Pos);
            triangleType = orderRemaingVerts(v2Pos, v3Pos, &vecPos);
        } else if(v2Pos[1] < v1Pos[1] && v2Pos[1] < v3Pos[1]){
            vecPos.push_back(v2Pos);
            triangleType = orderRemaingVerts(v1Pos, v3Pos, &vecPos);
        } else{
            vecPos.push_back(v3Pos);
            triangleType = orderRemaingVerts(v1Pos, v2Pos, &vecPos);
        }
        double leftLineGradient, rightLineGradient;
        double leftYIntercept, rightYIntercept;
        vector<int> tempVert;
        int v4x;
        int startX, endX;
        int xpos;
        
        calculateLineEquations(triangleType, &v4x, &tempVert,
            &leftLineGradient, &leftYIntercept, 
            &rightLineGradient, &rightYIntercept, 
            vecPos);

        xpos = vecPos[0][0];

        for(int y = vecPos[0][1]; y <= vecPos[2][1]; y++){
            if ((y - leftYIntercept) / leftLineGradient < (y - rightYIntercept) / rightLineGradient){
                startX = ((y - leftYIntercept) / leftLineGradient) + 0.5;
                endX = ((y - rightYIntercept) / rightLineGradient)+ 0.5;
            }else{
                startX = ((y - rightYIntercept) / rightLineGradient)+ 0.5;
                endX = ((y - leftYIntercept) / leftLineGradient)+ 0.5;
            }
            

            double edgeV1;
            double edgeV2;
            double edgeV3;

            for(xpos = startX; xpos<= endX; xpos++){
                edgeV1 = edgeFunction(v2, v3, xpos, y);
                edgeV2 = edgeFunction(v3, v1, xpos, y);
                edgeV3 = edgeFunction(v1, v2, xpos, y);
                if (edgeV1 >= 0 && edgeV2 >= 0 && edgeV3 >= 0){
                    if (xpos < 0 || y < 0){
                        continue;
                    }
                    if (xpos >= IMGWIDTH || y >= IMGHEIGHT){
                        continue;
                    }
                    double contribV1 = edgeV1 / area;
                    double contribV2 = edgeV2 / area;
                    double contribV3 = edgeV3 / area;

                    double faceZ = (contribV1 * vecPos[0][2]) + (contribV2 * vecPos[1][2]) + (contribV3 * vecPos[2][2]);
                    if (faceZ > (*framebuffer)[xpos+ 0.5][y].z && (*framebuffer)[xpos+ 0.5][y].z != 0.0){
                        continue;
                    }
                    vector<double> v1Col = v1->getSpecAndDiff();
                    vector<double> v2Col = v2->getSpecAndDiff();
                    vector<double> v3Col = v3->getSpecAndDiff();
                
                    (*framebuffer)[xpos][y].r = (contribV1 * v1Col[0]) + (contribV2 * v2Col[0]) + (contribV3 * v3Col[0]);
                    (*framebuffer)[xpos][y].g = (contribV1 * v1Col[1]) + (contribV2 * v2Col[1]) + (contribV3 * v3Col[1]);
                    (*framebuffer)[xpos][y].b = (contribV1 * v1Col[2]) + (contribV2 * v2Col[2]) + (contribV3 * v3Col[2]);
                    (*framebuffer)[xpos][y].z = faceZ;
                }
            }
            if (tempVert.size() > 0 && y == tempVert[1]){
                if (tempVert[0] > vecPos[1][0]){
                    //right side is longer
                    leftLineGradient = (((double)vecPos[1][1] - (double)vecPos[2][1]) / ((double)vecPos[1][0] - (double)vecPos[2][0]));
                    leftYIntercept = ((double)vecPos[1][1] - (leftLineGradient * (double)vecPos[1][0]));
                } else{
                    //left side is longer
                    rightLineGradient = (((double)vecPos[1][1] - (double)vecPos[2][1]) / ((double)vecPos[1][0] - (double)vecPos[2][0]));
                    rightYIntercept = ((double)vecPos[1][1] - (rightLineGradient * (double)vecPos[1][0]));
                }
            }
           
        }
        
        
    }
}

/*
 * Function: calculateLighting
 * ----------------------------
 *   Calculates all the lighting (the diffuse and the specular)
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 * 
 */
void calculateLighting(vector<vertex>* verticesPtr, vector<light>* lightsPtr, camera cam){
    calculateDiffuse(verticesPtr, lightsPtr);
    calculateSpecular(verticesPtr, lightsPtr, cam, 50.0);
}