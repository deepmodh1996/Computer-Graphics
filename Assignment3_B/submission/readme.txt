Roll Number : 140050002

Following advanced features are implemented :
	1. Refraction
	2. Elliptical Cylinder : One of the Quadric
	3. Triangles
	4. Constructive Solid Geometry
		Note : Only Primitives.[ch]pp are modified in case of CSG. Could not implement taking input due to time constraints. Approch used in taking Cylinder as input can be used here.

Bash Command to execute in data folder.
	5 sec for 		../trace  refraction.scd  ../images/refraction_image.png 2
	0.5 sec for 	../trace  cylinder.scd  ../images/cylinder_image.png 2
	5 sec for 		../trace  cylinder_refraction.scd  ../images/cylinder_refraction_image.png 2
	15 min for 		../trace  teapot.scd  ../images/teapot_image.png 2
	1.5 min for 	../trace  dunkit.scd  ../images/dunkit_image.png 2



------------------------------------------------------------------------------------------------------------------------------------

Refraction 

-------------------------------------------------------------------------------------------------------------------------------------

	Modified function :
		in main.cpp :
			RGB traceRay(Ray & ray, int depth)
		in types.hpp : in Class Ray
			double getN() const {return nVector_.back();} // vector of n ( refractive index)
			void pushN(double n);
		    void popN();

	Description : 

		in traceRay, if transparancy (mt) is non-zero, based on equations of refraction; calculated refracted ray and recursively call traceRay.
		To calculated refracted ray, we need refractive index of medium from which ray came or medium where ray will go after refraction.
		For that we need to store refractive indices. ( in nVector_ in class Ray)
		If refractive index of last entry in nVector_ is same as refractive index of primitive, then ray is going out of primitive, that is refractive index of outer medium is last entry in nVector_ after popping last entry( n of current medium).

	Images to check :
		refraction_image.png
		cylinder_refraction_image.png



-------------------------------------------------------------------------------------------------------------------------------------

Elliptical Cylinder : One of the Quadric

-------------------------------------------------------------------------------------------------------------------------------------

	Added function / class :
		SceneData.hpp
			class ParametricCylinder
				double getRadius(int time);
			    double getHeight(int time);
			    MaterialInfo getMaterial(int time);
			    ParametricCylinder()
			    ~ParametricCylinder()

		SceneGroup.[ch]pp
			bool SceneGroup::computeCylinder(double & radius, double & height, MaterialInfo & material, int time);

		SceneLoader.[ch]pp
			void setCylinderDefaults(SceneGroup * n);
			bool doCylinder(std::istream & str, std::string & name);

		Primitives.[ch]pp
			class Cylinder : public Primitive
				Cylinder()
				bool intersect(Ray & ray) const;
    			Vec3 calculateNormal(Vec3 const & position) const;

    Description :

    	Loading scene from .scd file to code :
    		Idea used is the same as that is used for loading sphere into code.
    		Scan through .scd file and store material information, radius, height and possible transformation, color and other properties of cylinder under map 'groups' as a parametricCylinder.
    		During importSceneToWorld function call in main.cpp, cylinder is created as an instance of Cylinder class which is subclass of Primitive. It is used to refer to cylinder in code.

    	Primitive :
    		In model Space, Cylinder is assumed to be of mentioned radius and height sitting on xy plane along positive z direction with centre of bottom circle as origin.
    		To calculate normal, check if point is on either of circle of cylinder or on curved surface.
    		To check if a ray intersects cylinder, consider atmost 4 intersections of ray with different part of cylinder.
    		One with each of the circle and two with curved surface.
    		For each possible intersection, check if such intersection is possible or not and if yes, if it has the minimum 't' parameter of ray.
    		To calculate point of intersection, qudratic and linear equations are solved.

    Images to check :
    	cylinder_image.png
    	cylinder_refraction_image.png



-------------------------------------------------------------------------------------------------------------------------------------

Triangles 

-------------------------------------------------------------------------------------------------------------------------------------
	
	Modified function :
		in Primitives.[ch]pp implemented functions of class Triangle :
			bool intersect(Ray & ray) const;
    		Vec3 calculateNormal(Vec3 const & position) const;

	Description :

		Same approch is used in case of Triangle as was used in case of Sphere.
		First calculate point of intersection using knowledge of vectors in 3-d cartesian coordinate system.
		Normal is constant for triangle, perpendicular to surface of triangle, outwards.
		To verify, teapot and dunkit are used which are given in skeleton code.

	Images to check :
		teapot_image.png
		dunkit_image.png



-------------------------------------------------------------------------------------------------------------------------------------

Constructive Solid Geometry

-------------------------------------------------------------------------------------------------------------------------------------
	
	Note : Only Primitives.[ch]pp are modified in case of CSG. Could not implement taking input due to time constraints. Approch used in taking Cylinder as input can be used here.

	Functions / class created :
		Primitives.[ch]pp :
			class CSG : public Primitive
				Primitive * A_; // first primitive ( can be another CSG Primitive)
			    Primitive * B_; // second primitive
			    int type_;
			    	type = 1 union
			    	type = 2 intersect
			    	type = 3 subtract : always A - B
				bool intersect(Ray & ray) const;
	    		Vec3 calculateNormal(Vec3 const & position) const;

	    	in all Primitives following functions are added. ( for 2d objects, CSG not considered hence returns false)
	    		bool contains(Vec3 const & position) const;
	    		bool lieOn(Vec3 const & position) const;

	Description :

		find intersect point of ray and CSG object :

			Let tA, tB be value of parameter 't' of ray, where ray intersects A and B respectively.
			Based on different cases, we determine intersection or recursively call this->intersect(newRay)
			for ease of explaining let step-0 be : check if minimum of tA, tB; say tC is less than minT of ray and return accordingly.
			in recursion, newRay has direction same as that of ray, but starts from 'tC' with slight offset.
			if ray intersects both objects;
				if type is union, do step-0.
				if type is intersection, point with parameter 'tC' on ray is contained in A and B, if yes; do step-0. else recurse.
				if type is subtraction, 
					if tA <= tB, do step-0
					else recurse.
			if intersects only A;
				if type is intersection, return false.
				else do step-0
			otherwise, if intersects only B;
				if type is union, do step-0
				else return false

		calculateNormal :
			position must lieOn surface of either A or B.
			return normal of A or B on which position lies.

		contains :
			based on type do boolean operation of contains(A), contains(B) by defination of type
				eg. in case of subraction, return ( (A->contains(position)) && (!(B->contains(position))) );

		lieOn :
			same as contains.
			based on type do boolean operation of lieOn(A), lieOn(B) by defination of type
				eg. in case of union, return ( (A->lieOn(position)) || (B->lieOn(position)) );
