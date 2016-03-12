#pragma once
#ifndef _OBJFRAME_
#define _OBJFRAME_

#include <vector>
#include <string>

class ObjFrame {
    
public:
    ObjFrame(int _res);
    void framerize(char *obj);
    int f2Ndx(float val);
    void calcFrameSize();
    
private:
    void writeFile(std::string filename);

    int getValAt(int x, int y, int z);
    void fillInFrame();
    void updateFrame(int x, int y, int z);

    void boxCenter(int x, int y, int z, float boxCntr[]);
   

    int res;
    float dx;
    std::vector<float> vertices;
    std::vector<int> faces;
    std::vector<int> frame;
};

#endif
