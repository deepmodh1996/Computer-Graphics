/*
 * Primitive.cpp
 *
 *  Created on: Feb 19, 2009
 *      Author: njoubert
 *      Modified: sidch
 */

#include "Primitives.hpp"

Primitive::Primitive(RGB const & c, Material const & m, Mat4 const & modelToWorld)
{
  c_ = c;
  m_ = m;
  modelToWorld_ = modelToWorld;
  worldToModel_ = modelToWorld.inverse();
}

Primitive::~Primitive()
{
}

Sphere::Sphere(double radius, RGB const & c, Material const & m, Mat4 const & modelToWorld): Primitive(c, m, modelToWorld)
{
  r_ = radius;
}

bool
Sphere::intersect(Ray & ray) const
{
  // TODO for 3a
  Ray rayModelSpace = Ray(ray);
  rayModelSpace.transform(this->worldToModel_);
  Vec3 e = rayModelSpace.start();
  Vec3 d = rayModelSpace.direction();
  double r = this->r_; 
  // center at origin
  double D = ( (d*e)*(d*e) ) - ( (d*d)*((e*e - r*r)) );// Discriminant
  if(D>=0.0){
    double t1 = ( (((-1.0)*d)*e) + sqrt(D) )/(d*d);
    double t2 = ( (((-1.0)*d)*e) - sqrt(D) )/(d*d);
    
    /* checking for specific conditions */
    double tMin = 0.0;
    double tMax = 0.0;
    if(t1>t2){tMin = t2;tMax = t1;}
    else{tMin = t1;tMax = t2;}

    /* modifying minT appropriately */
    if(tMin > 0.0){
      if(tMin < ray.minT()){
        ray.setMinT(tMin);
        return true;
      }else{return false;}
    }else{
      if(tMax>0.0){
         if(tMax < ray.minT()){
          ray.setMinT(tMax);
          return true;
        }else{return false;}
      }
    }

  }
  else{return false;}
  return false;
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

Vec3
Sphere::calculateNormal(Vec3 const & position) const
{
  // TODO for 3a
  return (worldToModel_.transpose() * (worldToModel_*position));
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

//=============================================================================================================================
// Triangle and other primitives are for Assignment 3b, after the midsem. Do not do this for 3a.
//=============================================================================================================================

Triangle::Triangle(Vec3 const & v0, Vec3 const & v1, Vec3 const & v2, RGB const & c, Material const & m,
                   Mat4 const & modelToWorld)
: Primitive(c, m, modelToWorld)
{
  verts[0] = v0;
  verts[1] = v1;
  verts[2] = v2;
}

bool
Triangle::intersect(Ray & ray) const
{
  // TODO for 3b, NOT 3a
  IMPLEMENT_ME(__FILE__, __LINE__);
}

Vec3
Triangle::calculateNormal(Vec3 const & position) const
{
  // TODO for 3b, NOT 3a
  IMPLEMENT_ME(__FILE__, __LINE__);
}
