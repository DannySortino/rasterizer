
#include "main.h"

using namespace std; 

/*
 * Function: calculateLightIntensity
 * ----------------------------
 *   Calculates the intensity of a light for a given fragment position
 * 
 *   vertexPt: A specific vertex
 *   lightPt: A specific light
 * 
 *   return: The calculated red, green and blue light intensity
 * 
 */
vector<double> calculateLightIntensity(vector<double> fragPos, vector<double> lightPosMat, vector<double> lightCol, double lightInt){
    double distance = getDistanceLength(fragPos, lightPosMat);
    double intensity = lightInt/pow(distance, 2.0);
    vector<double> RGB;
    RGB.push_back(lightCol[0] * intensity);
    RGB.push_back(lightCol[1] * intensity);
    RGB.push_back(lightCol[2] * intensity);

    return RGB;
}

/*
 * Function: addAmbient
 * ----------------------------
 *   Adds the ambient lighting to each fragment
 * 
 *   fragDiffCol: The diffuse colour of the fragment
 *   r: The red colour of the ambient
 *   g: The green colour of the ambient
 *   b: The blue colour of the ambient
 * 
 *   return: The calculated red, green and blue after adding ambient
 * 
 */
vector<double> addAmbient(vector<double> fragDiffCol, double r, double g, double b){
    vector<double> newColour;
    newColour.push_back(r * fragDiffCol[0]);
    newColour.push_back(g * fragDiffCol[1]);
    newColour.push_back(b * fragDiffCol[2]);
    return newColour;
}

/*
 * Function: calculateDiffuse
 * ----------------------------
 *   Calculates the diffuse for a fragment
 * 
 *   fragNormal: The normal of the fragment
 *   fragDiffCol: The diffuse colour of the fragment
 *   lightUnitVec: The light vector
 * 
 *   return: The calculated red, green and blue diffuse colour
 * 
 */
vector<double> calculateDiffuse(vector<double> fragNormal, vector<double> fragDiffCol, vector<double> lightUnitVec){
    vector<double> diff;
    double fragDot = dotProduct(fragNormal, lightUnitVec);
    if (fragDot < 0){
        fragDot = 0.0;
    }
    diff.push_back(fragDiffCol[0] * fragDot);
    diff.push_back(fragDiffCol[1] * fragDot);
    diff.push_back(fragDiffCol[2] * fragDot);
    return diff;
}

/*
 * Function: calculateSpecular
 * ----------------------------
 *   Calculates the specular for a fragment
 * 
 *   fragNormal: The normal of the fragment
 *   fragSpecCol: The specular colour of the fragment
 *   lightUnitVec: The light vector
 *   cameraPos: The position of the camera
 *   m: The specular exponent
 * 
 *   return: The calculated red, green and blue specular colour
 * 
 */
vector<double> calculateSpecular(vector<double> fragNormal, vector<double> fragSpecCol, vector<double> lightUnitVec, vector<double> cameraPos, double m){
    vector<double> halfwayVec = calculateHalfwayVector(lightUnitVec, cameraPos);
    double fragDot = dotProduct(fragNormal, halfwayVec);
    if (fragDot < 0){
        fragDot = 0.0;
    }
    double fragSpec = pow(fragDot, m);
    vector<double> spec;
    spec.push_back(fragSpecCol[0] * fragSpec);
    spec.push_back(fragSpecCol[1] * fragSpec);
    spec.push_back(fragSpecCol[2] * fragSpec);
        
    return spec;
}

/*
 * Function: calculateLighting
 * ----------------------------
 *   Calculates all the lighting (the diffuse and the specular)
 * 
 *   fragPos: The fragment position in world space
 *   fragNormal: The normal of the fragment
 *   fragDiffCol: The diffuse colour of the fragment
 *   fragSpecCol: The specular colour of the fragment
 *   fragLightVec: The light vector
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 * 
 *   return: The calculated red, green and blue for a fragment after adding the ambient, diffuse and specular
 * 
 */
vector<double> calculateLighting(vector<double> fragPos, vector<double> fragNormal, vector<double> fragDiffCol, vector<double> fragSpecCol, vector<double> fragLightVec, vector<light>* lightsPtr, camera cam){
    vector<double> colour = addAmbient(fragDiffCol, 0.3, 0.3, 0.3);

    for (int j = 0; j < lightsPtr->size(); j++){
        //Calculate light intensity
        vector<double> lightPosMat = (*lightsPtr)[j].getPosition();
        vector<double> lightCol = (*lightsPtr)[j].getRGB();
        double lightInt = (*lightsPtr)[j].getLightIntensity();
        vector<double> lightIntensity = calculateLightIntensity(fragPos, lightPosMat, lightCol, lightInt);

        //calculate diffuse
        vector<double> diff = calculateDiffuse(fragNormal, fragDiffCol, fragLightVec);

        //calculate specular
        vector<double> cameraPos = cam.getPosMat();
        vector<double> cameraUnitVec = getUnitLength(fragPos, cameraPos);
        vector<double> spec = calculateSpecular(fragNormal, fragSpecCol, fragLightVec, cameraUnitVec, 50.0);

        if (dotProduct(fragNormal, fragLightVec) > 0.0){
            colour[0] += (lightIntensity[0] * (diff[0] + spec[0]));
            colour[1] += (lightIntensity[1] * (diff[1] + spec[1]));
            colour[2] += (lightIntensity[2] * (diff[2] + spec[2]));
            if (colour[0] > 255.0){
                colour[0] = 255.0;
            }
            if (colour[1] > 255.0){
                colour[1] = 255.0;
            }
            if (colour[2] > 255.0){
                colour[2] = 255.0;
            }
        }            
    }
    return colour;
}

/*
 * Function: fragRasterize
 * ----------------------------
 *   Runs the main rasterising code, using a scanline approach, and calculaing the lighting for each fragment.
 * 
 *   verticesPtr: A pointer to a vector containing all vertices
 *   facesPtr: An array of the parameters passed to the application
 *   framebuffer: A 2D vector containing the frame buffer data that will be output to the image
 *   lightsPtr: A pointer to a vector of all lights
 *   cam: A class containing information about the camera
 * 
 */
void fragRasterize(vector<vertex>* verticesPtr, vector<face>* facesPtr, vector< vector<pixel> >* framebuffer, vector<light>* lightsPtr, camera cam){
    double near = 0.1;
    double far = 10.0;
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

                    //Calculate the position of the fragment in world space
                    vector<double> fragPos;
                    vector<double> v1WorldPos = v1->getPosMat();
                    vector<double> v2WorldPos = v2->getPosMat();
                    vector<double> v3WorldPos = v3->getPosMat();
                    fragPos.push_back((contribV1 * v1WorldPos[0]) + (contribV2 * v2WorldPos[0]) + (contribV3 * v3WorldPos[0]));
                    fragPos.push_back((contribV1 * v1WorldPos[1]) + (contribV2 * v2WorldPos[1]) + (contribV3 * v3WorldPos[1]));
                    fragPos.push_back((contribV1 * v1WorldPos[2]) + (contribV2 * v2WorldPos[2]) + (contribV3 * v3WorldPos[2]));

                    //Calculate the normal of the fragment
                    vector<double> fragNormal;
                    vector<double> v1Norm = v1->getNormal();
                    vector<double> v2Norm = v2->getNormal();
                    vector<double> v3Norm = v3->getNormal();
                    fragNormal.push_back((contribV1 * v1Norm[0]) + (contribV2 * v2Norm[0]) + (contribV3 * v3Norm[0]));
                    fragNormal.push_back((contribV1 * v1Norm[1]) + (contribV2 * v2Norm[1]) + (contribV3 * v3Norm[1]));
                    fragNormal.push_back((contribV1 * v1Norm[2]) + (contribV2 * v2Norm[2]) + (contribV3 * v3Norm[2]));
                    double distance = sqrt(pow(fragNormal[0], 2.0) + pow(fragNormal[1], 2.0) + pow(fragNormal[2], 2.0));
                    fragNormal[0] /= distance;
                    fragNormal[1] /= distance;
                    fragNormal[2] /= distance;

                    //Calculate the diffuse colour of the fragment
                    vector<double> fragDiffCol;
                    vector<double> v1DiffCol = v1->getColour();
                    vector<double> v2DiffCol = v2->getColour();
                    vector<double> v3DiffCol = v3->getColour();
                    fragDiffCol.push_back((contribV1 * v1DiffCol[0]) + (contribV2 * v2DiffCol[0]) + (contribV3 * v3DiffCol[0]));
                    fragDiffCol.push_back((contribV1 * v1DiffCol[1]) + (contribV2 * v2DiffCol[1]) + (contribV3 * v3DiffCol[1]));
                    fragDiffCol.push_back((contribV1 * v1DiffCol[2]) + (contribV2 * v2DiffCol[2]) + (contribV3 * v3DiffCol[2]));

                    //Calculate the specular colour of the fragment
                    vector<double> fragSpecCol;
                    vector<double> v1SpecCol = v1->getSpecColour();
                    vector<double> v2SpecCol = v2->getSpecColour();
                    vector<double> v3SpecCol = v3->getSpecColour();
                    fragSpecCol.push_back((contribV1 * v1SpecCol[0]) + (contribV2 * v2SpecCol[0]) + (contribV3 * v3SpecCol[0]));
                    fragSpecCol.push_back((contribV1 * v1SpecCol[1]) + (contribV2 * v2SpecCol[1]) + (contribV3 * v3SpecCol[1]));
                    fragSpecCol.push_back((contribV1 * v1SpecCol[2]) + (contribV2 * v2SpecCol[2]) + (contribV3 * v3SpecCol[2]));

                    //Calculate the fragment light vector
                    vector<double> fragLightVec;
                    vector<double> lightPosMat = (*lightsPtr)[0].getPosition();
                    vector<double> v1fragLight = getUnitLength(v1WorldPos, lightPosMat);
                    vector<double> v2fragLight = getUnitLength(v2WorldPos, lightPosMat);
                    vector<double> v3fragLight = getUnitLength(v3WorldPos, lightPosMat);
                    fragLightVec.push_back((contribV1 * v1fragLight[0]) + (contribV2 * v2fragLight[0]) + (contribV3 * v3fragLight[0]));
                    fragLightVec.push_back((contribV1 * v1fragLight[1]) + (contribV2 * v2fragLight[1]) + (contribV3 * v3fragLight[1]));
                    fragLightVec.push_back((contribV1 * v1fragLight[2]) + (contribV2 * v2fragLight[2]) + (contribV3 * v3fragLight[2]));

                    distance = sqrt(pow(fragLightVec[0], 2.0) + pow(fragLightVec[1], 2.0) + pow(fragLightVec[2], 2.0));
                    fragLightVec[0] /= distance;
                    fragLightVec[1] /= distance;
                    fragLightVec[2] /= distance;


                    vector<double> colour = calculateLighting(fragPos, fragNormal, fragDiffCol, fragSpecCol, fragLightVec, lightsPtr, cam);

                    (*framebuffer)[xpos][y].r = round(colour[0]);
                    (*framebuffer)[xpos][y].g = round(colour[1]);
                    (*framebuffer)[xpos][y].b = round(colour[2]);
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