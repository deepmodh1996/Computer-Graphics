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
}

Vec3
Sphere::calculateNormal(Vec3 const & position) const
{
  return (worldToModel_.transpose() * (worldToModel_*position));
}

bool 
Sphere::contains(Vec3 const & position) const
{
  return (worldToModel_*position).length() <= r_;
}

bool 
Sphere::lieOn(Vec3 const & position) const
{
  Vec3 positionModel = worldToModel_*position;
  double x = positionModel.x();
  double y = positionModel.y();
  double z = positionModel.z();
  if( (x*x + y*y + z*z) == (r_*r_) ){return true;}
  return false;
}

//=============================================================================================================================
// Triangles
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
  Ray rayModelSpace = Ray(ray);
  rayModelSpace.transform(this->worldToModel_);

  // get triangle edge vectors and plane normal
  Vec3 u = verts[1] - verts[0];
  Vec3 v = verts[2] - verts[0];
  Vec3 n = u ^ v;              // cross product

  Vec3 dir = rayModelSpace.direction();              // ray direction vector
  Vec3 w = rayModelSpace.start() - verts[0];
  double a = (-1.0)*(n*w);
  double b = n*dir;
  if (b == 0.0) {return false;}     // ray is  parallel to triangle plane

  // get intersect point of ray with triangle plane
  double r = a / b;
  if (r < 0.0){return false;}                    // ray goes away from triangle
  //     return 0;                   // => no intersect
  // // for a segment, also test if (r > 1.0) => no intersect

  Vec3 pt = rayModelSpace.start() + r * dir;            // intersect point of ray and plane

  // is pt inside T?
  double uu = u*u;
  double uv = u*v;
  double vv = v*v;
  // w change
  Vec3 x = pt - verts[0];
  double xu = x*u;
  double xv = x*v;
  double D = (uv * uv) - (uu * vv);

  // get and test parametric coords
  double s, t;
  s = (uv * xv - vv * xu) / D;
  if (s < 0.0 || s > 1.0)         // I is outside T
      return false;
  t = (uv * xu - uu * xv) / D;
  if (t < 0.0 || (s + t) > 1.0)  // I is outside T
      return false;

  double tMin = ray.minT();
  if(r < tMin){
    ray.setMinT(r);
    return true;  
  }
  return false;
}

Vec3
Triangle::calculateNormal(Vec3 const & position) const
{
  Vec3 u = verts[1] - verts[0];
  Vec3 v = verts[2] - verts[0];
  Vec3 n = u ^ v;
  return (worldToModel_.transpose() *n);
}

bool 
Triangle::contains(Vec3 const & position) const
{
  return false;
}

bool 
Triangle::lieOn(Vec3 const & position) const
{
  return false;
}


//////////////////////////////////////////////////////////////////////////
///       CSG
////////////////////////////////////////////////////////////////////////////

CSG::CSG(Primitive * A, Primitive * B, int type, RGB const & c, Material const & m, Mat4 const & modelToWorld): Primitive(c, m, modelToWorld)
{
  A_ = A;
  B_ = B;
  type_ = type;
}

bool
CSG::intersect(Ray & ray) const
{
  Ray rayModelSpaceA = Ray(ray);
  rayModelSpaceA.transform(this->worldToModel_);
  Ray rayModelSpaceB = Ray(ray);
  rayModelSpaceB.transform(this->worldToModel_);

  bool isIntersectA = A_->intersect(rayModelSpaceA);
  bool isIntersectB = B_->intersect(rayModelSpaceB);

  if(isIntersectA && isIntersectB){
    double tA = rayModelSpaceA.minT();
    double tB = rayModelSpaceB.minT();
    double tMin = 0.0;
    if(tA <= tB){tMin = tA;}
    else{tMin = tB;}

    if(type_ == 1){
      if(tMin < ray.minT()){
        ray.setMinT(tMin);
        return true;
      }else{return false;}
    }
    else{
      if(type_ == 2){
        Vec3 pos = ray.start() + tMin * ray.direction();
        // pos in argument of contains must be in world coordinates
        if( (A_->contains(pos)) && (B_->contains(pos))){
          if(tMin < ray.minT()){
            ray.setMinT(tMin);
            return true;
          }else{return false;}
        }
        else{
           Ray rayTemp = Ray(rayModelSpaceA);
           Vec3 posTB = rayModelSpaceA.start() + tMin * rayModelSpaceA.direction();
           Ray newRay = rayTemp.fromOriginAndDirection( posTB + (0.001 * rayModelSpaceA.direction()), rayModelSpaceA.direction() );
           return this->intersect(newRay);
         }
      }
    else{ // type == 3
      if(tMin == tA){// tA <= tB
        if(tMin < ray.minT()){
          ray.setMinT(tMin);
          return true;
        }else{return false;}
      }
      else{
        Ray rayTemp = Ray(rayModelSpaceA);
        Vec3 posTB = rayModelSpaceA.start() + tMin * rayModelSpaceA.direction();
        Ray newRay = rayTemp.fromOriginAndDirection( posTB + (0.001 * rayModelSpaceA.direction()), rayModelSpaceA.direction() );
        return this->intersect(newRay);
      }
    }

    return false;
    }
  }

  if(isIntersectA){
    if(type_ == 2){return false;}
    else{
      if(rayModelSpaceA.minT() < ray.minT()){
        ray.setMinT(rayModelSpaceA.minT());
        return true;
      }else{return false;}
    }

    return false;
  }

  if( (type_ == 1) && ( rayModelSpaceB.minT() < ray.minT() ) ){
    ray.setMinT(rayModelSpaceB.minT());
    return true;
  }
    
   return false;


}

// position lies on CSG thus on surface of either A or B; using this fact
 Vec3
 CSG::calculateNormal(Vec3 const & position) const
 {
    if(A_->lieOn(position)){return A_->calculateNormal(position);}
    return B_->calculateNormal(position);
}


 bool 
 CSG::contains(Vec3 const & position) const
 {
  if(type_ == 1){
    return ( (A_->contains(position)) || (B_->contains(position)) );
  }
  else{
    if(type_ == 2){
      return ( (A_->contains(position)) && (B_->contains(position)));
    }
    else{
      return ( (A_->contains(position)) && (!(B_->contains(position))) );
    }

  }
  return false;
  }

   bool 
 CSG::lieOn(Vec3 const & position) const
 {
  if(type_ == 1){
    return ( (A_->lieOn(position)) || (B_->lieOn(position)) );
  }
  else{
    if(type_ == 2){
      return ( (A_->lieOn(position)) && (B_->lieOn(position)));
    }
    else{
      return ( (A_->lieOn(position)) && (!(B_->lieOn(position))) );
    }

  }
  return false;
  }






/////////////////////////////////////////////////////////////////////////
///////       Cylinder
/////////////////////////////////////////////////////////////////////////

Cylinder::Cylinder(double radius, double height, RGB const & c, Material const & m, Mat4 const & modelToWorld): Primitive(c, m, modelToWorld)
{
  r_ = radius;
  h_ = height;
}

bool
Cylinder::intersect(Ray & ray) const
{
  // TODO for 3a
  Ray rayModelSpace = Ray(ray);
  rayModelSpace.transform(this->worldToModel_);
  Vec3 e = rayModelSpace.start();
  Vec3 d = rayModelSpace.direction();

  double t1 = -1.0;
  double t2 = -1.0;
  if(d.z()!=0){
    t1 = ( h_ - e.z() ) / d.z();
    t2 = ( 0.0 - e.z() ) / d.z();
    Vec3 I1 = (e + t1*d); // Intersection 1
    Vec3 I2 = (e + t2*d); // Intersection 2
    if((I1.x()*I1.x() + I1.y()*I1.y()) > r_*r_){t1 = -1.0;}
    if((I2.x()*I2.x() + I2.y()*I2.y()) > r_*r_){t2 = -1.0;}
  }

  double ex = e.x();
  double ey = e.y();
  double dx = d.x();
  double dy = d.y();
  double a = ((dx*dx)+(dy*dy));
  double b = 2.0*((ex*dx)+(ey*dy));
  double c = ((ex*ex)+(ey*ey)-(r_*r_));
  double D = (b*b) - (4.0*a*c);

  double t3 = -1.0;
  double t4 = -1.0;

  if(D==0.0){
    t3 = ((-1.0)*b) / (2*a);
    double zTemp = (e + t3*d).z();
    if((zTemp<0.0) || (zTemp>h_)){t3 = -1.0;}
  }

  if(D>0.0){
    t3 = (((-1.0)*b)-sqrt(D)) / (2*a);
    t4 = (((-1.0)*b)+sqrt(D)) / (2*a);
    double zTemp = (e + t3*d).z();
    if((zTemp<0.0) || (zTemp>h_)){t3 = -1.0;}
    zTemp = (e + t4*d).z();
    if((zTemp<0.0) || (zTemp>h_)){t4 = -1.0;}
  }

  double tMin = ray.minT();
  if((t1<tMin) && (t1>=0.0)){tMin = t1;}
  if((t2<tMin) && (t2>=0.0)){tMin = t2;}
  if((t3<tMin) && (t3>=0.0)){tMin = t3;}
  if((t4<tMin) && (t4>=0.0)){tMin = t4;}

  if(tMin<ray.minT()){ray.setMinT(tMin);return true;}

  return false; 
}

Vec3
Cylinder::calculateNormal(Vec3 const & position) const
{
  Vec3 positionModel = worldToModel_*position;
  Vec3 normalModel = Vec3(0.0,0.0,0.0);
  double x = positionModel.x();
  double y = positionModel.y();
  double z = positionModel.z();
  double epsilon = 0.000001;
  if((z<h_+epsilon)&&(z>h_-epsilon)){
    normalModel = Vec3(0.0,0.0,1.0);
  }
  else if((z<epsilon)&&(z>(-1.0*epsilon))){
    normalModel = Vec3(0.0,0.0,-1.0);
  }
  else{
    normalModel = Vec3(x,y,0.0);
  }
  return (worldToModel_.transpose() * (normalModel));
  
}

bool 
Cylinder::contains(Vec3 const & position) const
{
  Vec3 positionModel = worldToModel_*position;
  double x = positionModel.x();
  double y = positionModel.y();
  double z = positionModel.z();
  if( (x*x + y*y) > (r_*r_)){return false;}
  if( (z>h_) || (z<0.0) ){return false;}
  return true;
}

bool 
Cylinder::lieOn(Vec3 const & position) const
{
  Vec3 positionModel = worldToModel_*position;
  double x = positionModel.x();
  double y = positionModel.y();
  double z = positionModel.z();
  if( ((x*x + y*y) == (r_*r_)) && ((z>h_) || (z<0.0)) ){return true;}
  return false;
}