#include "ObjFrame.h"
#include "TriBoxOverlap.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <cmath>

using namespace std;

ObjFrame::ObjFrame(int _res) :
    res(_res)
{
}

void ObjFrame::framerize(char *obj) 
{

    frame.resize(res*res*res, 0);

    // Read Sequence
    ifstream in;
    in.open(obj);
    if (!in.good()) {
        cout << "Cannot read " << obj << endl;
        return;
    }
    
    string line;
    string id;

    float x, y, z;
    int f1, f2, f3;
    
    while(1) {
        getline(in, line);
        
        if (in.eof()) {
            break;
        } 

        // Skip empty lines
		if(line.size() < 2) {
			continue;
		}
		// Skip comments
		if(line.at(0) == '#') {
			continue;
		}

        stringstream ss(line);
        
        ss >> id;
        
        if (id == "v") {
            
            ss >> x >> y >> z;
            
            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
            
        } else if (id == "f") {
            
            ss >> f1 >> f2 >> f3;

            faces.push_back(f1);
            faces.push_back(f2);
            faces.push_back(f3);
            
        } else {
            continue;
        }
    }

    in.close();

    calcFrameSize();
    fillInFrame();

    string str = string(obj);
    string newF = str.substr(str.find("/")+1, str.find("."));
    newF = newF.substr(0, newF.find("."));

    writeFile("frames/" + newF + ".frame");
}

void ObjFrame::writeFile(string filename)
{
    
    const char *f = filename.c_str();
    ofstream out(f);
    if(!out.good()) {
        cout << "Could not open " << filename << endl;
        return;
    }
    
    for (int k = 0; k < res; k++) {
        for (int j = 0; j < res; j++) {
            for (int i = 0; i < res; i++) {
                out << getValAt(i, j, k) << " ";
            }
            out << endl;
        }
        out << endl;
    }
    
    out.close();
}

void ObjFrame::calcFrameSize() 
{
    float best = 0.0;

    for (int i = 0; i < vertices.size(); i++) {
        best = max(vertices[i], abs(best));
    }

    dx = (2 * best) / res;
}

int ObjFrame::getValAt(int x, int y, int z)
{
    return frame[x + y*res + z*res*res];
}
   

void ObjFrame::fillInFrame()
{

    float boxHlfSize[3] = {dx/2, dx/2, dx/2};
    float triverts[3][3];
    int vert;
    float bxCntr[3];
    
    // Conduct check on each face
    for (int i = 0; i < faces.size()/3; i++) {
        // Acquire the face vertices
        for (int j = 0; j < 3; j++) {
            // starting index of vert
            vert = faces[i*3 + j];
            for (int k = 0; k < 3; k++) {
                triverts[j][k] = vertices[vert + k];
            }
        }       
        
        // calc bounding box of cubes to check intersections of face with
        float maxX = max(max(triverts[0][0], triverts[1][0]), triverts[2][0]);
        float minX = min(min(triverts[0][0], triverts[1][0]), triverts[2][0]);
        float maxY = max(max(triverts[0][1], triverts[1][1]), triverts[2][1]);
        float minY = min(min(triverts[0][1], triverts[1][1]), triverts[2][1]);
        float maxZ = max(max(triverts[0][2], triverts[1][2]), triverts[2][2]);
        float minZ = min(min(triverts[0][2], triverts[1][2]), triverts[2][2]);
        
        int bmaxX = f2Ndx(maxX);
        int bminX = f2Ndx(minX);
        int bmaxY = f2Ndx(maxY);
        int bminY = f2Ndx(minY);
        int bmaxZ = f2Ndx(maxZ);
        int bminZ = f2Ndx(minZ);

        //cout << bmaxX << " " << bminX << " " << bmaxY << " " << bminY << " " << bmaxZ << " " << bminZ << endl;

        // Check whether triangle intersects the cubes
        for (int x = bminX; x < bmaxX; x++) {
            for (int y = bminY; y < bmaxY; y++) {
                for (int z = bminZ; z < bmaxZ; z++) {
                    boxCenter(x, y, z, bxCntr);
                    
                    if (triBoxOverlap(bxCntr, boxHlfSize, triverts) == 1) { 

                        cout << "X: " << x << " Y: " << y << " Z: " << x << endl;
                        // Fill frame in 
                        updateFrame(x, y, z);
                    }
                }
            }
        }
    }
}

         

void ObjFrame::updateFrame(int x, int y, int z)
{
    frame[x + y*res + z*res*res] = 1;
}

int ObjFrame::f2Ndx(float val)
{
    // Shift vertex component since obj centered at origin
    val += dx * (res/2);
    
    // Make sure calculated ndx is within bounds of frame
    int ndx = (val > dx * res) ? (res - 1) : val / dx; 

    return ndx; 
}

void ObjFrame::boxCenter(int x, int y, int z, float boxCntr[]) {
    float boxHalf = dx / 2;
    
    // Calcs pos in 1st quadrant then shifts to center
    boxCntr[0] = (x*dx + dx/2) - 1.0;
    boxCntr[1] = (y*dx + dx/2) - 1.0;
    boxCntr[2] = (z*dx + dx/2) - 1.0;
}
