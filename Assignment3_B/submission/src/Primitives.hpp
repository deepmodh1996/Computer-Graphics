/*
 * Primitive.hpp
 *
 *  Created on: Feb 19, 2009
 *      Author: njoubert
 *      Modified: sidch
 */

#ifndef __Primitive_hpp__
#define __Primitive_hpp__

#include "Globals.hpp"

/** Interface for a scene primitive (e.g. a sphere). */
class Primitive
{
  protected:
    Primitive(RGB const & c, Material const & m, Mat4 const & modelToWorld);  // primitives are an abstract class

  public:
    /** Destructor. */
    virtual ~Primitive();

    /**
     * Checks for intersection with the given ray. If there is a valid intersection which is smaller than the ray's current
     * minimum hit time (ray.minT()), then updates ray.minT() to the new hit time (as a multiple of the ray length). The ray is
     * specified in world space.
     *
     * @return True if there was a valid intersection AND the minimum hit time was lowered.
     *
     * !!! REMEMBER TO TAKE THE PRIMITIVE TRANSFORM (modelToWorld_/worldToModel_) INTO ACCOUNT !!!
     */
    virtual bool intersect(Ray & ray) const = 0;

    /**
     * Calculates the normal for the given position on this primitive. You may assume the position is actually on the
     * primitive. The position is specified in world space.
     *
     * !!! REMEMBER TO TAKE THE PRIMITIVE TRANSFORM (modelToWorld_/worldToModel_) INTO ACCOUNT !!!
     */
    virtual Vec3 calculateNormal(Vec3 const & position) const = 0;

    // for 3d objects position in world coordinates is contained within object
    // false for 2d objects vacuously
    virtual bool contains(Vec3 const & position) const = 0;

    /** Set the primitive's color. */
    void setColor(RGB const & c) { c_ = c; }

    /** Set the primitive's material. */
    void setMaterial(Material const & m) { m_ = m; }

    /** Get the primitive's color. */
    RGB const & getColor() const { return c_; }

    /** Get the primitive's material. */
    Material const & getMaterial() const { return m_; }

    virtual bool lieOn(Vec3 const & position) const = 0;

  protected:
    Mat4 modelToWorld_;
    Mat4 worldToModel_;

  private:
    RGB c_;
    Material m_;
};

/** A sphere primitive. */
class Sphere : public Primitive
{
  public:
    /** Constructor. */
    Sphere(double radius, RGB const & c, Material const & m, Mat4 const & modelToWorld);

    bool intersect(Ray & ray) const;
    Vec3 calculateNormal(Vec3 const & position) const;
    bool contains(Vec3 const & position) const;
    bool lieOn(Vec3 const & position) const;

  private:
    double r_;
};

/** A Cylinder primitive. */
class Cylinder : public Primitive
{
  public:
    /** Constructor. */
    Cylinder(double radius, double height, RGB const & c, Material const & m, Mat4 const & modelToWorld);

    bool intersect(Ray & ray) const;
    Vec3 calculateNormal(Vec3 const & position) const;
    bool contains(Vec3 const & position) const;
    bool lieOn(Vec3 const & position) const;

  private:
    double r_;
    double h_;
};

//=============================================================================================================================
// Triangle and other primitives are for Assignment 3b, after the midsem. Do not do this for 3a. The skeleton is included here
// as a preview.
//=============================================================================================================================

/** A triangle primitive. */
class Triangle : public Primitive
{
  public:
    /** Constructor. */
    Triangle(Vec3 const & v0, Vec3 const & v1, Vec3 const & v2, RGB const & c, Material const & m, Mat4 const & modelToWorld);

    bool intersect(Ray & ray) const;
    Vec3 calculateNormal(Vec3 const & position) const;
    bool contains(Vec3 const & position) const;
    bool lieOn(Vec3 const & position) const;

  private:
    Vec3 verts[3];
};

/** A CSG primitive. */
class CSG : public Primitive
{
  public:
    /** Constructor. */
    // Discription of A and B in world Space is taken as Their discription in modelSpace of CSG,
    // Thus modelToWorld of CSG transform world discription of A and B according to transformation of CSG object
    // in generel set modelToWorld to Identity
    
    // material, color is initialized to all values as zero (Actual material and color will be of A, B according to type)

    
    // type = 1 union
    // type = 2 intersect
    // type = 3 subtract : always A - B
    CSG(Primitive * A, Primitive * B, int type, RGB const & c, Material const & m, Mat4 const & modelToWorld);
    bool intersect(Ray & ray) const;
    Vec3 calculateNormal(Vec3 const & position) const;
    bool contains(Vec3 const & position) const;
    bool lieOn(Vec3 const & position) const;

  private:
    Primitive * A_;
    Primitive * B_;
    int type_;
};


#endif  // __Primitive_hpp__
