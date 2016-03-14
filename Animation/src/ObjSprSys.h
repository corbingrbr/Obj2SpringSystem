#pragma once
#ifndef _OBJSPRSYS_
#define _OBJSPRSYS_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

class Particle;
class Spring;
class MatrixStack;
class Program;
class Shape;

class ObjSprSys
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
        ObjSprSys(std::string fileName, Eigen::Vector3d begPos, double scale, double mass, 
                  double stiffness, const Eigen::Vector2d &damping,
                  std::shared_ptr<Shape> sphere);
	virtual ~ObjSprSys();
	
    void genParticlesNSprings();
	void init();

	void tare();
	void reset();
	void step(double h, const Eigen::Vector3d &grav);
	
	void drawParticles(std::shared_ptr<MatrixStack> MV, std::shared_ptr<Program> p) const;
    void drawSprings();

    Eigen::Vector3d getPos();
    double getRadius();
    void colliding(std::shared_ptr<Particle> p);

private:
    
    int gRef(int x, int y, int z);
    int gNdx(int x, int y, int z, int dim);
    void makeSpring(int ndx, int x, int y, int z, int dim, 
                    int os, int ondx, int ofx, int ofy, int ofz, 
                    std::vector<std::shared_ptr<Particle> > &frame, 
                    std::vector<bool> &s);

    void genSprings(int x, int y, int z, int dim, std::vector<std::shared_ptr<Particle> > &frame, std::vector<bool> &s);
    


	int n;
    double scale;
    double mass;
    double stiffness;
    Eigen::Vector3d pos;
    std::string fileName;
    Eigen::Vector3d begPos;
	Eigen::Vector2d damping;
    std::shared_ptr<Shape> sphere;
	std::vector< std::shared_ptr<Particle> > particles;
	std::vector< std::shared_ptr<Spring> > springs;
    
    const int NS = 26; // Potential number of springs a particle can have
    
    // Reference for setting up springs
    const int ref[3][3][3] = {{ {0,1,2},{3,4,5},{6,7,8} },
                              { {9,10,11},{12,-1,13},{14,15,16} },
                              { {17,18,19},{20,21,22},{23,24,25} }};
	
	Eigen::VectorXd v;
	Eigen::VectorXd f;
	Eigen::MatrixXd M;
	Eigen::MatrixXd K;
};

#endif
