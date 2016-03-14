#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "ObjSprSys.h"
#include "Shape.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

extern bool drawSprings;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{
	// Units: meters, kilograms, seconds
	h = 1e-2;
	
	grav << 0.0, -4.0, 0.0;
	
	double mass = 0.01;
	double stiffness = 1e2;
    Vector2d damping(0.0, 1.0);

    auto sphere = make_shared<Particle>(sphereShape);

	sphereShape = make_shared<Shape>();
    sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
    
    // Load Cube 
    Vector3d begPos(0,1.0,0);
    objSprSys = make_shared<ObjSprSys>(RESOURCE_DIR + "cube.frame", begPos, mass, stiffness, damping, sphereShape);

}

void Scene::init()
{
	sphereShape->init();
	objSprSys->init();
}

void Scene::tare()
{
	objSprSys->tare();
}

void Scene::reset()
{
	t = 0.0;
	objSprSys->reset();
}

void Scene::step()
{
	t += h;

	// Simulate the cloth
	objSprSys->step(h, grav); // Possibly include shapes in scene for collisions*/
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if (!drawSprings) {
        glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
    }
	objSprSys->draw(MV, prog);
}
