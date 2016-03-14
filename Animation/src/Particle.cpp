#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <stdlib.h>

#include "Particle.h"
#include "Shape.h"
#include "Program.h"
#include "MatrixStack.h"

using namespace std;
using namespace Eigen;

Particle::Particle() :
	r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
	fixed(true)
{
	
}

Particle::Particle(const shared_ptr<Shape> s) :
    r(1.0),
	m(1.0),
	i(-1),
	x(0.0, 0.0, 0.0),
	v(0.0, 0.0, 0.0),
    fixed(true)
{
    
}

Particle::Particle(const shared_ptr<Shape> s, double r, Vector3d x0) :
	m(.05),
	i(-1),
	v(0.0, 0.0, 0.0),
	fixed(false),
	sphere(s)
{
    this->r = r;
    this->x0 = x0;
    this->x = x0;
    this->v << 0.0, 0.0, 0.0;
}


Particle::~Particle()
{
}

void Particle::tare()
{
	x0 = x;
	v0 = v;
}

void Particle::reset()
{
	x = x0;
	v << 0.0, 0.0, 0.0;
}

void Particle::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog) const
{
	if(sphere) {

		MV->pushMatrix();
		MV->translate(Eigen::Vector3f(x(0), x(1), x(2)));
		MV->scale(r);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, MV->topMatrix().data());
		sphere->draw(prog);
		MV->popMatrix();
	}
}

void Particle::colliding(shared_ptr<Particle> p) 
{
    if (!fixed && (x - p->x).norm() <= r + p->r) {
        
        Vector3d v1, v2, v1x, v1y, v2x, v2y, dx; 
        float m1, m2, d1, d2;
        
        dx = x - p->x;
        dx.normalize();
        
        //x = dx * p->r;

        v1 = v;
        
        d1 = dx.dot(v1);
        v1x = dx * d1;
        v1y = v1 - v1x;
        m1 = m;
        
        dx = dx*-1;
        v2 = p->v;
        d2 = dx.dot(v2);
        v2x = dx * d2;
        v2y = v2 - v2x;
        m2 = p->m;
        
        v = v1x*(m1-m2)/(m1+m2) + v2x*(2*m2)/(m1+m2) + v1y;    
        
    }
}
