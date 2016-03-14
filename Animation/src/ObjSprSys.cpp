#include "ObjSprSys.h"
#include "Particle.h"
#include "Spring.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"
#include "GLSL.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

shared_ptr<Spring> createSpring(const shared_ptr<Particle> p0, const shared_ptr<Particle> p1, double E)
{
    auto s = make_shared<Spring>(p0, p1);
    s->E = E;
    Vector3d x0 = p0->x;
    Vector3d x1 = p1->x;
    Vector3d dx = x1 - x0;
	s->L = dx.norm();
	return s;
}

ObjSprSys::ObjSprSys(string fileName, Vector3d begPos, double scale, double mass, 
                     double stiffness, const Vector2d &damping,
                     shared_ptr<Shape> sphere)
{
	assert(mass > 0.0);
	assert(stiffness > 0.0);
	
    this->fileName = fileName;
    this->pos = begPos;
    this->begPos = begPos;
	this->damping = damping;
    this->scale = scale;
	this->mass = mass;
    this->stiffness = stiffness;	
    this->sphere = sphere;
}

ObjSprSys::~ObjSprSys()
{
}

void ObjSprSys::genParticlesNSprings()
{
    // Read through file
    ifstream in;
    in.open(fileName);
    if (!in.good()) {
        cout << "Cannot read " << fileName << endl;
        return;
    }
   
    string line;   
    int dim, val;
   
    n = 0;
    
    in >> dim; // Dimension of frame x y z   
    getline(in, line); // Discard rest of line
    
    double os = -0.5*scale; // origin shift
    double d = 1.0 / dim * scale; // particle diameter
    double r = d / 2; // particle radius
    double m = mass / (dim * dim * dim); // Particle mass    
                                                vector<shared_ptr<Particle> > frame; // 3D vector for holding particle structure
    frame.assign(dim*dim*dim, NULL); // Initialize all to null pointers 

    for (int k = 0; k < dim; k++) {
        for (int j = 0; j < dim; j++) {
            for (int i = 0; i < dim; i++) {
                in >> val;
                if (val == 1) {
                    Vector3d x;
                    x << (i+1)*d-r + os, (j+1)*d-r + os, (k+1)*d-r + os;
                    x += begPos;
                    auto p = make_shared<Particle>(sphere, r, x);
                    p->m = m;
                    p->i = n;
                    n += 3;
                    particles.push_back(p);
                    frame[i + j*dim + k*dim*dim] = p;
                } 
            }
            getline(in, line); // Discard
        }
        getline(in, line); // Discard
    }

    in.close();
    
    vector<bool> s; // For keeping track of spring connections between particles
    s.assign(dim*dim*dim*NS, false);

    for (int z = 0; z < dim; z++) {
        for (int y = 0; y < dim; y++) {
            for (int x = 0; x < dim; x++) {
                if (frame[gNdx(x, y, z, dim)] != NULL) { // Find first particle in structure
                    genSprings(x, y, z, dim, frame, s);          // genSprings from that point recursively
                    return;
                }
            }
        }
    }
}

void ObjSprSys::init()
{
    // constructs particles and springs from frame file
    genParticlesNSprings();

	// Build system matrices and vectors given number of particles
	M.resize(n,n);
	K.resize(n,n);
	v.resize(n);
	f.resize(n);
}

void ObjSprSys::tare()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->tare();
	}
}

void ObjSprSys::reset()
{
	for(int k = 0; k < (int)particles.size(); ++k) {
		particles[k]->reset();
	} 
    v.setZero();
}

void ObjSprSys::step(double h, const Vector3d &grav)
{
    // Update So it uses sparse matrices
    
	M.setZero();
	K.setZero();
	v.setZero();
	f.setZero();

    MatrixXd I(3,3);
    I << MatrixXd::Identity(3,3);

    int ndx;
    
    for (unsigned int i = 0; i < particles.size(); i++) {
        // Add particle data to memory
        ndx = particles[i]->i;        
        M.block<3,3>(ndx,ndx) = particles[i]->m * I;
        v.segment<3>(ndx) << particles[i]->v;
        f.segment<3>(ndx) << particles[i]->m * grav;
    }
    
    // Add forces from springs to particles. 
    for (unsigned int i = 0; i < springs.size(); i++) {
        Vector3d dx = springs[i]->p1->x - springs[i]->p0->x;

        double l = dx.norm(); // magnitude        
        double E = springs[i]->E;
        double L = springs[i]->L;

        // Indices of the particles that compose the spring 
        int pndx0 = springs[i]->p0->i;
        int pndx1 = springs[i]->p1->i;
        
        double exp = (l - L) / l;
        
        MatrixXd Ks(3,3);
        
        double dxtdx = dx.transpose() * dx;
        Ks = E/(l*l) * ((1-exp)*(dx * dx.transpose()) + exp*dxtdx*I);
        
        K.block<3,3>(pndx0, pndx0) += Ks;  // Top Left, Bottom Right
        K.block<3,3>(pndx1, pndx1) += Ks; 
        
        K.block<3,3>(pndx0, pndx1) += -Ks; // Top Right, Bottom Left
        K.block<3,3>(pndx1, pndx0) += -Ks; 
        
        
        dx.normalize(); // normalized vector for direction of force
        Vector3d fs = E*(l-L)*dx;
       
        f.segment<3>(pndx0) += fs;
        f.segment<3>(pndx1) += -fs; 
    }
    

    MatrixXd A(n,n);
    VectorXd b(n);
    

    // Implicit Euler: M + h^2 * K v^(k+1) = M*v^k + h*f^k
    // A = M + h^2 * K
    MatrixXd D(n,n);
    D = damping(0)*h*M + damping(1)*h*h*K;
    A = M + D;
    b = M*v + h*f;
    
    VectorXd s(n);

    // Solve Ax = b 
    s = A.ldlt().solve(b);

    // Update velocities
    for (unsigned int i = 0; i < particles.size(); i++) {
        ndx = particles[i]->i;
        particles[i]->v = s.segment<3>(ndx);
        particles[i]->x += h*particles[i]->v;
    }
}


void ObjSprSys::drawParticles(shared_ptr<MatrixStack> MV, const shared_ptr<Program> p) const
{
        // Draw particles
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->draw(MV, p);
    }
}

void ObjSprSys::drawSprings()
{
    // Draw spring system
    for (unsigned int i = 0; i < springs.size(); i++) {
        springs[i]->draw();
    }
}

int ObjSprSys::gRef(int x, int y, int z) 
{
    return ref[z+1][y+1][x+1];
}

int ObjSprSys::gNdx(int x, int y, int z, int dim)
{
    return x + y*dim + z*dim*dim;
}


void ObjSprSys::makeSpring(int ndx, int x, int y, int z, int dim, int os, int ondx, int ofx, int ofy, int ofz, vector<shared_ptr<Particle> > &frame, vector<bool> &s)
{
     s[os+gRef(ofx, ofy, ofz)] = true;
     s[ondx*NS + gRef(-ofx, -ofy, -ofz)] = true;
     springs.push_back(createSpring(frame[ndx], frame[ondx], stiffness));
     genSprings(x+ofx, y+ofy, z+ofz, dim, frame, s);    
}

void ObjSprSys::colliding(shared_ptr<Particle> object)
{
    // Run through particles
    for (unsigned int i = 0; i < particles.size(); i++) {
        particles[i]->colliding(object);
    }
}

void ObjSprSys::genSprings(int x, int y, int z, int dim, vector<shared_ptr<Particle> > &frame, vector<bool> &s) 
{ 
    // Frame to check if particle exists in a location
    // Springs to check if a spring exists between 

    bool zgz = z > 0;
    bool zld = z < dim - 1;
    bool ygz = y > 0;
    bool yld = y < dim - 1;
    bool xgz = x > 0;
    bool xld = x < dim - 1;

    // Check if search is possible in a direction
    // Check if particle exists
    // Check if spring exists between the two
    // If not then mark springs between the two () -> createSpring

    assert(x >= 0 && x < dim && y >= 0 && y < dim && z >= 0 && z < dim);

    
    int ndx = gNdx(x,y,z, dim); // index of particle
    int ondx; // index of other particle
    int os = ndx * NS; // offset into springs
    
    // Chxeck z plane behind
    if (zgz) {
        // Check y plane below
        if (ygz) {    

            ondx = gNdx(x-1, y-1, z-1, dim);

            if (xgz && frame[ondx] && !s[os+gRef(-1, -1, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, -1, -1, frame, s); }
            
            ondx = gNdx(x, y-1, z-1, dim);
            if (frame[ondx] && !s[os+gRef(0, -1, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, -1, -1, frame, s); }
            
            ondx = gNdx(x+1, y-1, z-1, dim);
            if (xld && frame[ondx] && !s[os+gRef(1, -1, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, -1, -1, frame, s); }
        }
        
        //Check same y plane
        ondx = gNdx(x-1, y, z-1, dim);
        if (xgz && frame[ondx] && !s[os+gRef(-1, 0, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, 0, -1, frame, s); }
        
        ondx = gNdx(x, y, z-1, dim);
        if (frame[ondx] && !s[os+gRef(0, 0, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, 0, -1, frame, s); }
        
        ondx = gNdx(x+1, y, z-1, dim);
        if (xld && frame[ondx] && !s[os+gRef(1, 0, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 0, -1, frame, s); }
        
        //Check y plane above
        if (yld) {
            ondx = gNdx(x-1, y+1, z-1, dim);
            if (xgz && frame[ondx] && s[os+gRef(-1, 1, -1)]) { 
                makeSpring(ndx, x, y, z, dim, os, ondx, -1, 1, -1, frame, s); }
            
            ondx = gNdx(x, y+1, z-1, dim);
            if (frame[ondx] && s[os+gRef(0, 1, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, 1, -1, frame, s); }
                
            ondx = gNdx(x+1, y+1, z-1, dim);
            if (xld && frame[ondx] && s[os+gRef(1, 1, -1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 1, -1, frame, s); }
        }
    }
  
    //Check same z plane
    
    // Check y plane below
        if (ygz) {
            ondx = gNdx(x-1, y-1, z, dim);
            if (xgz && frame[ondx] && !s[os+gRef(-1, -1, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, -1, 0, frame, s); }
                
            ondx = gNdx(x, y-1, z, dim);
            if (frame[ondx] && !s[os+gRef(0, -1, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, -1, 0, frame, s); }
            
            ondx = gNdx(x+1,y-1,z, dim);
            if (xld && frame[ondx] && !s[os+gRef(1, -1, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, -1, 0, frame, s); }
        }
        
        //Check same y plane
        ondx = gNdx(x-1, y, z, dim);
        if (xgz && frame[ondx] && !s[os+gRef(-1, 0, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, 0, 0, frame, s); }
            
        ondx = gNdx(x+1, y, z, dim);
        if (xld && frame[ondx] && !s[os+gRef(1, 0, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 0, 0, frame, s); }

        
        //Check y plane above
        if (yld) {
            ondx = gNdx(x-1, y+1, z, dim);
            if (xgz && frame[ondx] && !s[os+gRef(-1, 1, 0)]) {
                makeSpring(ndx, x, y, z, dim, os, ondx, -1, 1, 0, frame, s); }
   
            ondx = gNdx(x, y+1, z, dim);
                if (frame[ondx] && !s[os+gRef(0, 1, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, 1, 0, frame, s); }
            
            ondx = gNdx(x+1, y+1, z, dim);
            if (xld && frame[ondx] && !s[os+gRef(1, 1, 0)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 1, 0, frame, s); }
        }
  
    
    // Check z plane in front
    if (zld) {
        // Check y plane below
        if (ygz) {
            ondx = gNdx(x-1, y-1, z+1, dim);
            if (xgz && frame[ondx] && !s[os+gRef(-1, -1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, -1, 1, frame, s); }

            ondx = gNdx(x, y-1, z+1, dim);
            if (frame[ondx] && !s[os+gRef(0, -1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, -1, 1, frame, s); }
            
            ondx = gNdx(x+1, y-1, z+1, dim);
            if (xld && frame[ondx] && !s[os+gRef(1, -1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, -1, 1, frame, s); }
        }
        
        //Check same y plane
        ondx = gNdx(x-1,y,z+1, dim);
        if (xgz && frame[ondx] && !s[os+gRef(-1, 0, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, 0, 1, frame, s); }
  
        ondx = gNdx(x,y,z+1, dim);
        if (frame[ondx] && !s[os+gRef(0, 0, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, 0, 1, frame, s); }
       
        ondx = gNdx(x+1,y, z+1, dim);
        if (xld && frame[ondx] && !s[os+gRef(1, 0, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 0, 1, frame, s); }
        
        //Check y plane above
        if (yld) {
            ondx = gNdx(x-1, y+1 , z+1, dim);
            if (xgz && frame[ondx] && !s[os+gRef(-1, 1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, -1, 1, 1, frame, s); }
            
            ondx = gNdx(x, y+1,z+1, dim);
            if (frame[ondx] && !s[os+gRef(0, 1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 0, 1, 1, frame, s); }
                
            ondx = gNdx(x+1, y+1, z+1, dim);
            if (xld && frame[ondx] && !s[os+gRef(1, 1, 1)]) { makeSpring(ndx, x, y, z, dim, os, ondx, 1, 1, 1, frame, s); }
        }     
    }
}
