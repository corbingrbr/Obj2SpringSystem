#include <iostream>

#include "Scene.h"
#include "Particle.h"
#include "ObjSprSys.h"
#include "Shape.h"
#include "Program.h"

using namespace std;
using namespace Eigen;

extern bool drawSprngs;

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
	h = 2e-3;
	
	grav << 0.0, -9.8, 0.0;
	
	double mass = 0.01;
	//double stiffness = 1.2e2;
	double stiffness = 5.0; // 10.0
    Vector2d damping(0.0, 1.0);
    double scale = 0.1;

    auto sphere = make_shared<Particle>(sphereShape);

	sphereShape = make_shared<Shape>();
    sphereShape->loadMesh(RESOURCE_DIR + "sphere2.obj");
    
    // Load Cube 
    Vector3d begPos(.15,1.0,0);
    objSprSys = make_shared<ObjSprSys>(RESOURCE_DIR + "cube.frame", begPos, scale, mass, stiffness, damping, sphereShape);

    // Create Sphere object for collisions
    double r = .15;
    begPos << 0, r, 0; 
    sphrr = make_shared<Particle>(sphereShape, r, begPos);
    sphrr->fixed = true;

    r = .2;
    begPos << .5, -.5, 0; 
    sphrr2 = make_shared<Particle>(sphereShape, r, begPos);
    sphrr2->fixed = true;

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
	if (!drawSprngs) {
        objSprSys->drawParticles(MV, prog);
    } 

    glUniform3fv(prog->getUniform("kdFront"), 1, Vector3f(1.0, 1.0, 1.0).data());
    sphrr->draw(MV, prog);
    sphrr2->draw(MV, prog);
}

void Scene::drawSprings(shared_ptr<MatrixStack>MV, const shared_ptr<Program> prog) const
{
    objSprSys->drawSprings();
}

void Scene::collisionDetection() 
{
    if (isColliding(sphrr, objSprSys)) {
        objSprSys->colliding(sphrr);
    }
    
    if (isColliding(sphrr2, objSprSys)) {
        objSprSys->colliding(sphrr2);
    }  
}

bool Scene::isColliding(shared_ptr<Particle> sphere, shared_ptr<ObjSprSys> obj)
{
    return true;
}
