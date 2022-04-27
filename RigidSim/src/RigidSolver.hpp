// ----------------------------------------------------------------------------
// RigidSolver.hpp
//
//  Created on: 18 Dec 2020
//      Author: Kiwon Um
//        Mail: kiwon.um@telecom-paris.fr
//
// Description: Simple Rigid Body Solver (DO NOT DISTRIBUTE!)
//
// Copyright 2020 Kiwon Um
//
// The copyright to the computer program(s) herein is the property of Kiwon Um,
// Telecom Paris, France. The program(s) may be used and/or copied only with
// the written permission of Kiwon Um or in accordance with the terms and
// conditions stipulated in the agreement/contract under which the program(s)
// have been supplied.
// ----------------------------------------------------------------------------

#ifndef _RIGIDSOLVER_HPP_
#define _RIGIDSOLVER_HPP_

#include <glm/ext/matrix_transform.hpp>

#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include "Quaternion.hpp"

struct BodyAttributes {
  BodyAttributes() :
    X(0, 0, 0), R(Mat3f::I()), P(0, 0, 0), L(0, 0, 0),
    V(0, 0, 0), omega(0, 0, 0), F(0, 0, 0), tau(0, 0, 0), Q(Quaternion(1.0, 0.0, 0.0, 0.0)) {}

  glm::mat4 worldMat() const
  {
    return glm::mat4(           // column-major
      R(0,0), R(1,0), R(2,0), 0,
      R(0,1), R(1,1), R(2,1), 0,
      R(0,2), R(1,2), R(2,2), 0,
      X[0],   X[1],   X[2],   1);
  }

  Real M;                       // mass
  Mat3f I0, I0inv;              // inertia tensor and its inverse in body space

  // rigid body state
  Vec3f X;                      // position
  Mat3f R;                      // rotation
  Vec3f P;                      // linear momentum
  Vec3f L;                      // angular momentum

  // auxiliary quantities
  Mat3f Iinv;                   // inverse of inertia tensor
  Vec3f V;                      // linear velocity
  Vec3f omega;                  // angular velocity
  Quaternion Q;

  // force and torque
  Vec3f F;                      // force
  Vec3f tau;                    // torque

  // mesh's vertices in body space
  std::vector<Vec3f> vdata0;
};

class Box : public BodyAttributes {
public:
  explicit Box(
    const Real w=1.0, const Real h=1.0, const Real d=1.0, const Real dens=10.0,
    const Vec3f v0=Vec3f(0, 0, 0), const Vec3f omega0=Vec3f(0, 0, 0)) :
    width(w), height(h), depth(d)
  {
    V = v0;                     // initial velocity
    omega = omega0;             // initial angular velocity

    M = w * h * d * dens;
    I0 = Mat3f(
      M*(h*h+d*d)/12, 0, 0,
      0, M*(w*w+d*d)/12, 0,
      0, 0, M*(w*w+h*h)/12
    );
    I0inv = I0.inverse();

    // vertices data (8 vertices)
    vdata0.push_back(Vec3f(-0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h, -0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h, -0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h, -0.5*d));

    vdata0.push_back(Vec3f(-0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w, -0.5*h,  0.5*d));
    vdata0.push_back(Vec3f( 0.5*w,  0.5*h,  0.5*d));
    vdata0.push_back(Vec3f(-0.5*w,  0.5*h,  0.5*d));
  }

  // rigid body property
  Real width, height, depth;
};

class RigidSolver {
public:
  explicit RigidSolver(
    BodyAttributes *body0=nullptr, const Vec3f g=Vec3f(0, 0, 0), const Real drag=0.0) :
    body(body0), _g(g), _step(0), _sim_t(0), _drag(drag) {}

  void init(BodyAttributes *body0)
  {
    body = body0;
    _step = 0;
    _sim_t = 0;
  }

  void step(const Real dt)
  {
    std::cout << "t=" << _sim_t << " (dt=" << dt << ")" << std::endl;

    computeForceAndTorque();

    // linear momentum
    body->P += dt * body->F;
    body->V = body->P / body->M;
    body->X += dt * body->V;

    // rotation
    body->L += dt * body->tau;
    body->Iinv = body->R * body->I0inv * body->R.transposed();
    body->omega = body->Iinv * body->L;
    body->Q += Quaternion(0, body->omega.x, body->omega.y, body->omega.z) * body->Q * dt / 2;
    body->R = body->Q.normalize().ToMat();


    ++_step;
    _sim_t += dt;
  }

  BodyAttributes *body;

private:
  void computeForceAndTorque()
  {
    // reset force and torque
    body->tau = Vec3f(0, 0, 0);
    body->F = Vec3f(0, 0, 0);

    // add gravity
    body->F += _g * body->M;

    // add drag force and torque
    for(Vec3f vdata: body->vdata0){
      // compute force on each vertices with body speed + vertex local speed
      Mat3f localSpace = crossProductMatrix(body->R * vdata); // precalculate local space vector ; and its cross product matrix
      Vec3f force = _drag * (body->V  - localSpace * body->omega) / body->vdata0.size() / 2; // force applied to the vertex
      body->F -= force; // substract force
      body->tau -= localSpace * force; // substract it to torque
    }

    // instant forces and torque at step 1
    if(_step==1) {
        body->F = Vec3f(0.15f, 0.25f, 0.03f);

        body->tau = crossProductMatrix(body->R * body->vdata0[0]) * Vec3f(0.15f, 0.25f, 0.03f);
    }
  }

  // simulation parameters
  Real _drag;
  Vec3f _g;                     // gravity
  tIndex _step;                 // step count
  Real _sim_t;                  // simulation time
};

#endif  /* _RIGIDSOLVER_HPP_ */
