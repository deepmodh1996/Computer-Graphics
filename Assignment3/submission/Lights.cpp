/*
 * Lights.cpp
 *
 *  Created on: Feb 19, 2009
 *      Author: njoubert
 *      Modified: sidch
 */

#include "Lights.hpp"

Light::Light()
{
  RGB black(0, 0, 0);
  illumination_ = black;
  falloff_ = 0;
  angular_falloff_ = 0;
  dead_distance_ = 1;
}

Light::Light(RGB const & illumination)
{
  illumination_ = illumination;
}

Light::Light(RGB const & illumination, double falloff, double dead_distance)
{
  illumination_ = illumination;
  falloff_ = falloff;
  dead_distance_ = dead_distance;
}

Light::~Light()
{}

RGB Light::getColor() const
{
  return illumination_;
}

RGB Light::getColor(Vec3 const & p) const
{
  return illumination_;
}

void
Light::setColor(RGB const & c)
{
  illumination_ = c;
}

AmbientLight::AmbientLight()
{
  // intentionally empty
}

AmbientLight::AmbientLight(RGB const & illumination) : Light(illumination)
{
  // intentionally empty
}

Vec3
AmbientLight::getIncidenceVector(Vec3 const & position) const
{
  throw "AMBIENT LIGHTS DO NOT HAVE A SENSE OF DIRECTION OR POSITION`";
}

Ray AmbientLight::getShadowRay(Vec3 const & position, bool & use_dist) const
{
  throw "AMBIENT LIGHTS DO NOT HAVE A SENSE OF DIRECTION OR POSITION";
}

PointLight::PointLight(RGB const & illumination) : Light(illumination)
{
  // intentionally empty
}

PointLight::PointLight(RGB const & illumination, double falloff, double dead_distance)
: Light(illumination, falloff, dead_distance)
{
  // intentionally empty
}

RGB
PointLight::getColor(Vec3 const & p) const
{
  // TODO for 3a
  double dist = Vec3(p - (this->pos_)).length();
    double scale = 1.0/(pow(dist + (this->dead_distance_), this->falloff_));
    return (this->illumination_)*scale;
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

void
PointLight::setPosition(Vec3 const & pos)
{
  pos_ = pos;
}

Vec3
PointLight::getIncidenceVector(Vec3 const & position) const
{
  // TODO for 3ag
  Vec3 returnVec = (this->pos_) - position;
  return returnVec.normalize();
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

Ray
PointLight::getShadowRay(Vec3 const & position, bool & use_dist) const
{
  // TODO for 3a
  // IMPLEMENT_ME(__FILE__, __LINE__);
  use_dist = true;
  Ray returnRay;
  return returnRay.fromOriginAndDirection(position, pos_ - position);
}

DirectionalLight::DirectionalLight(RGB const & illumination) : Light(illumination)
{
  // intentionally empty
}

void
DirectionalLight::setDirection(Vec3 const & dir)
{
  dir_ = dir;
  dir_.normalize();
}

Vec3
DirectionalLight::getIncidenceVector(Vec3 const & position) const
{
  // TODO for 3a
  return (-1.0)*(this->dir_);
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

Ray
DirectionalLight::getShadowRay(Vec3 const & position, bool & use_dist) const
{
  // TODO for 3a
  // IMPLEMENT_ME(__FILE__, __LINE__);
  use_dist = false;
  Ray returnRay;
  return returnRay.fromOriginAndDirection(position, (-1.0)*(this->dir_));
}
