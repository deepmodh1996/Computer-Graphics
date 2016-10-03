/*
 * World.cpp
 *
 *  Created on: Feb 19, 2009
 *      Author: njoubert
 *      Modified: sidch
 */

#include "World.hpp"

World::World()
{
}

World::~World()
{
  // TODO Auto-generated destructor stub
}

Primitive *
World::intersect(Ray & r) const
{
	Primitive *intersectPrimitivePtr = NULL;
	long long vSize = primitives_.size();
	for(long long i=0;i<vSize;i++){
		if((primitives_[i])->intersect(r)){intersectPrimitivePtr = primitives_[i];}
	}
	return intersectPrimitivePtr;
  // IMPLEMENT_ME(__FILE__, __LINE__);
}

void
World::addPrimitive(Primitive * p)
{
  primitives_.push_back(p);
}

void
World::addLight(Light * l)
{
  lights_.push_back(l);
}

void
World::setAmbientLightColor(RGB ambientColor)
{
  ambientLight_.setColor(ambientColor);
}

RGB
World::getAmbientLightColor() const
{
  return ambientLight_.getColor();
}

void
World::printStats() const
{
  std::cout << "World data:" << std::endl;
  std::cout << " primitives: " << primitives_.size() << std::endl;
  std::cout << " lights: " << lights_.size() << std::endl;
}
