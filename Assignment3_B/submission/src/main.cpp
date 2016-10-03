#include "Globals.hpp"
#include "View.hpp"
#include "World.hpp"
#include "Frame.hpp"
#include "Lights.hpp"
#include "core/Scene.hpp"

using namespace std;

//****************************************************
// Global Variables
//****************************************************

Scene * scene = NULL;
World * world = NULL;
View * view = NULL;
Mat4 viewToWorld = identity3D();
Frame * frame = NULL;
int max_trace_depth = 2;

// Get the shaded appearance of the primitive at a given position, as seen along a ray. The returned value should be the sum of
// the shaded colors w.r.t. each light in the scene. DO NOT include the result of recursive raytracing in this function, just
// use the ambient-diffuse-specular formula. DO include testing for shadows, individually for each light.
RGB
getShadedColor(Primitive const & primitive, Vec3 const & pos, Ray const & ray)
{
  // IMPLEMENT_ME(__FILE__, __LINE__);

  // Use the "world" global variable to access the lights in the input file.
  //  Please be sure to read the following classes in Types.hpp:
  //   - RGB
  //   - Material
  //   - Ray

  RGB C = primitive.getColor(); // color of object
  RGB S =  ( primitive.getMaterial().getMSM() * C ) + ( (1-primitive.getMaterial().getMSM()) * RGB(1.0,1.0,1.0) ); // specular highlight
  RGB ambientColor =  primitive.getMaterial().getMA() * (  C * world->getAmbientLightColor());
  RGB finalColor = ambientColor;

  std::vector<Light *>::const_iterator it = world->lightsBegin();
  while(it!=world->lightsEnd()){

    /* getting material property */
    double ml = primitive.getMaterial().getML();
    double ms = primitive.getMaterial().getMS();
    double msp = primitive.getMaterial().getMSP();
    
    RGB I = (*it)->getColor(pos);// intensity of infalling light
    
    /* getting all unit vectors */
    Vec3 ICap = (*it)->getIncidenceVector(pos);
    Vec3 tempVec = ICap.normalize();
    
    Vec3 uCap = ray.direction();
    tempVec = uCap.normalize();
    
    Vec3 nCap = primitive.calculateNormal(pos);
    tempVec = nCap.normalize();
    
    Vec3 rCap = (-1.0 * ICap) + ( (2*(ICap*nCap)) * nCap );
    double dotproduct = ((-1.0)*rCap)*uCap; 
    double exponent = pow(max(dotproduct,0.0) , msp );
    
    /* checking for shadow */
    bool use_dist = false;
    Ray shadowRay = (*it)->getShadowRay(pos + (0.001 * nCap),use_dist);
    Primitive * shadowPrimitiveIntersect = world->intersect(shadowRay);
    if(shadowPrimitiveIntersect!=NULL){
      
      if(!use_dist){
        if(shadowRay.minT() != std::numeric_limits<double>::infinity()){std::advance (it,1);continue;}
      }
      else{
        if( (shadowRay.minT() <= 1.0) && (shadowRay.minT()>0.0) ){std::advance (it,1);continue;}
      }
    }

    /* calculating color */
    finalColor = finalColor +   ( (ml * (C * I)) * max(ICap*nCap, 0.0) ) + ( (ms * (S * I)) * exponent );
  
    std::advance (it,1);
  }
  return finalColor;
}

// Raytrace a single ray backwards into the scene, calculating the total color (summed up over all reflections/refractions) seen
// along this ray.
RGB
traceRay(Ray & ray, int depth)
{
  if (depth > max_trace_depth){
    return RGB(0.0, 0.0, 0.0);
  }

  // IMPLEMENT_ME(__FILE__, __LINE__);

  // Use the "world" global variable to access the primitives in the input file.
  //  IMPORTANT:
  //  Please start all bounce rays at a small non-zero t value such as 0.001 - this has the effect of slightly offsetting
  //  bounce rays from the surface they're bouncing from, and prevents bounce rays from being occluded by their own surface.

  Primitive *intersectPrimitivePtr = world->intersect(ray);

  if(intersectPrimitivePtr == NULL){return RGB(0.0, 0.0, 0.0);}
  else{
    
    Vec3 pos = ray.getPos(ray.minT());
    RGB directColor =  getShadedColor(*intersectPrimitivePtr, pos, ray);
    
    /* calculating bouncing ray */
    Vec3 nCap = intersectPrimitivePtr->calculateNormal(pos);
    Vec3 vtemp = nCap.normalize();
    Vec3 incidentDirection = ray.direction();
    vtemp = incidentDirection.normalize();
    
    /* calculating bounce color */
    Vec3 bounceDirection =  incidentDirection -  (  (2*(incidentDirection*nCap)) * nCap);
    Ray rayTemp;
    Ray bounceRay = rayTemp.fromOriginAndDirection( pos + (0.001 * nCap), bounceDirection);  
    RGB bounceColor = traceRay(bounceRay, depth + 1);

    RGB refractedColor = RGB(0.0,0.0,0.0);
    double mt = intersectPrimitivePtr->getMaterial().getMT();
    if(mt!=0){  
      double n2 = intersectPrimitivePtr->getMaterial().getMTN();
      double n1 = ray.getN();
      if(n1 == n2){
        ray.popN();
        n2 = ray.getN();
      }
      /* calculating refracting ray*/
      double cos1 = (-1.0)*(incidentDirection*nCap);
      double sin2Sqr = (1.0-(cos1*cos1)) * ( (n1*n1 / (n2*n2) ));
      if(sin2Sqr<=1.0){
        double cos2 = sqrt(1.0 - sin2Sqr);  
        Vec3 refractDirection = (n1/n2)*(incidentDirection) + ( ((n1/n2)*cos1)- cos2)*nCap;
        vtemp = refractDirection.normalize();
        Ray refractedRay = rayTemp.fromOriginAndDirection( pos + (0.001 * (-1.0) * nCap), refractDirection );
        refractedRay.pushN(n2);
        refractedColor = mt*traceRay(refractedRay, depth+1);
      }
    }

    RGB finalBounceColor = (intersectPrimitivePtr->getMaterial().getMR() * intersectPrimitivePtr->getColor() * bounceColor);
    return ( directColor + finalBounceColor + refractedColor);
    
  }
}

// Main rendering loop.
void
renderWithRaytracing()
{
  Sample sample;   // Point on the view being sampled.
  Ray ray;         // Ray being traced from the eye through the point.
  RGB c;           // Color being accumulated per pixel.

  int const rpp = view->raysPerPixel();

  for (int yi = 0; yi < view->height(); ++yi)
  {
    for (int xi = 0; xi < view->width(); ++xi)
    {
      c = RGB(0, 0, 0);
      for (int ri = 0; ri < rpp; ++ri)
      {
        view->getSample(xi, yi, ri, sample);
        ray = view->createViewingRay(sample);  // convert the 2d sample position to a 3d ray
        ray.transform(viewToWorld);            // transform this to world space
        c += traceRay(ray, 0);
      }

      frame->setColor(sample, c / (double)rpp);
    }
  }
}

// This traverses the loaded scene file and builds a list of primitives, lights and the view object. See World.hpp.
void
importSceneToWorld(SceneInstance * inst, Mat4 localToWorld, int time)
{
  if (inst == NULL)
    return;

  Mat4 nodeXform;
  inst->computeTransform(nodeXform, time);
  localToWorld = localToWorld * nodeXform;
  SceneGroup * g = inst->getChild();

  if (g == NULL)   // for example if the whole scene fails to load
  {
    std::cout << "ERROR: We arrived at an instance with no child?!" << std::endl;
    return;
  }

  int ccount = g->getChildCount();

  for (int i = 0; i < ccount; i++)
  {
    importSceneToWorld(g->getChild(i), localToWorld, 0);
  }

  CameraInfo f;

  if (g->computeCamera(f, time))
  {
    viewToWorld = localToWorld;

    if (view != NULL)
      delete view;

    Vec3 eye(0.0, 0.0, 0.0);
    Vec3 LL(f.sides[FRUS_LEFT], f.sides[FRUS_BOTTOM], -f.sides[FRUS_NEAR]);
    Vec3 UL(f.sides[FRUS_LEFT], f.sides[FRUS_TOP], -f.sides[FRUS_NEAR]);
    Vec3 LR(f.sides[FRUS_RIGHT], f.sides[FRUS_BOTTOM], -f.sides[FRUS_NEAR]);
    Vec3 UR(f.sides[FRUS_RIGHT], f.sides[FRUS_TOP], -f.sides[FRUS_NEAR]);
    view = new View(eye, LL, UL, LR, UR, IMAGE_WIDTH, IMAGE_HEIGHT, RAYS_PER_PIXEL_EDGE);
  }

  LightInfo l;

  if (g->computeLight(l, time))
  {
    if (l.type == LIGHT_AMBIENT)
    {
      RGB amb = world->getAmbientLightColor();
      world->setAmbientLightColor(amb + l.color);
    }
    else if (l.type == LIGHT_DIRECTIONAL)
    {
      DirectionalLight * li = new DirectionalLight(l.color);
      Vec3 dir(0, 0, -1);
      li->setDirection(localToWorld * dir);
      world->addLight(li);
    }
    else if (l.type == LIGHT_POINT)
    {
      PointLight * li = new PointLight(l.color, l.falloff, l.deadDistance);
      Vec3 pos(0, 0, 0);
      li->setPosition(localToWorld * pos);
      world->addLight(li);
    }
    else if (l.type == LIGHT_SPOT)
    {
      throw "oh no";
    }
  }

  double r;
  MaterialInfo m;

  if (g->computeSphere(r, m, time))
  {
    Material mat(m.k[0], m.k[1], m.k[2], m.k[3], m.k[4], m.k[MAT_MS], m.k[5], m.k[6]);
    Sphere * sph = new Sphere(r, m.color, mat, localToWorld);
    world->addPrimitive(sph);
  }

  double h;
  if (g->computeCylinder(r, h, m, time))
  {
    Material mat(m.k[0], m.k[1], m.k[2], m.k[3], m.k[4], m.k[MAT_MS], m.k[5], m.k[6]);
    Cylinder * cyl = new Cylinder(r, h, m.color, mat, localToWorld);
    world->addPrimitive(cyl);
  }

  TriangleMesh * t;

  if (g->computeMesh(t, m, time))
  {
    Material mat(m.k[0], m.k[1], m.k[2], m.k[3], m.k[4], m.k[MAT_MS], m.k[5], m.k[6]);

    for (vector<MeshTriangle *>::iterator it = t->triangles.begin(); it != t->triangles.end(); ++it)
    {
      Triangle * tri = new Triangle(
        t->vertices[ (**it).ind[0] ]->pos,
        t->vertices[ (**it).ind[1] ]->pos,
        t->vertices[ (**it).ind[2] ]->pos,
        m.color, mat, localToWorld);
      world->addPrimitive(tri);
    }
  }

  std::cout << "Imported scene file" << std::endl;
}

int
main(int argc, char ** argv)
{
  if (argc < 3)
  {
    std::cout << "Usage: " << argv[0] << " scene.scd output.png [max_trace_depth]" << std::endl;
    return -1;
  }

  if (argc >= 4)
    max_trace_depth = atoi(argv[3]);

  cout << "Max trace depth = " << max_trace_depth << endl;

  // Load the scene from the disk file
  scene = new Scene(argv[1]);

  // Setup the world object, containing the data from the scene
  world = new World();
  importSceneToWorld(scene->getRoot(), identity3D(), 0);
  world->printStats();

  // Set up the output framebuffer
  frame = new Frame(IMAGE_WIDTH, IMAGE_HEIGHT);

  // Render the world
  renderWithRaytracing();

  // Save the output to an image file
  frame->save(argv[2]);
  std::cout << "Image saved!" << std::endl;
  std::cout<< max(-0.232, 0.0)<<std::endl;
}
