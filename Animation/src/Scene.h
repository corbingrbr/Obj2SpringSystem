#pragma once
#ifndef __Scene__
#define __Scene__

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class ObjSprSys;
class Particle;
class MatrixStack;
class Program;
class Shape;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void tare();
	void reset();
	void step();
    void collisionDetection();
	
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog) const;
	void drawSprings(std::shared_ptr<MatrixStack>MV, const std::shared_ptr<Program> prog) const;
	double getTime() const { return t; }
	
private:
    
    bool isColliding(std::shared_ptr<Particle> sphere, std::shared_ptr<ObjSprSys> obj);

	double t;
	double h;
	Eigen::Vector3d grav;
    std::shared_ptr<Particle> sphrr;
    std::shared_ptr<Particle> sphrr2;
	
	std::shared_ptr<Shape> sphereShape;
    std::shared_ptr<ObjSprSys> objSprSys;

};

#endif
