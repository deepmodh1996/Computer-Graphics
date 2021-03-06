#include "Algebra3.hpp"
#include "Image.hpp"
#include "LineSegment.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/*****************************************************************************
Morphing consists of distorting each image, and blending the two results. I
suggest you first write a blend function that generates each output pixel as a
linear interpolation of the two input pixels at the same location. Once you
have this working, you can blend the two input images and see whether things
are working correctly. This should produce an image that looks like both input
images are semi-transparent and overlaid on one another.

Once you have blend working, you can attempt the much harder distortion
function. This takes a list of line segments at the start and end of the time,
and an image corresponding to the image at the start. It then interpolates
each segment as a combination of the start and end segment. This new list of
segments is used to distort the image according to the algorithm described in
Feature-Based Image Metamorphosis, Beier & Neely, 1992.

You will use the distort function on both input images to create a distorted
version of both images at the same point in time, which means that img1 needs
to be moved from 0 to t (total amount t), while img2 needs to be moved from 1
to t (total amount 1 - t). Since we want to write a single function to do
this, you need to swap the segments and the t value you pass to the distort
function for img2 versus img1.

A good way to debug the distort function is to return the result of the
distort on img1 or img2 without blending the two, so you can see if the
resulting distort works. Good luck!
*****************************************************************************/


/**
 * Use bilinear interpolation to get the color at an image position \a loc with real-valued coordinates, by interpolating from
 * the 4 surrounding pixels, and return the result in \a sampled_color. In the diagram below, P denotes image pixel positions
 * (integer row and column), and L is a position to be sampled. In this case the color should be interpolated from the corner
 * pixels of the top-left square.
 *
 * NOTES:
 * - Do NOT use half-integer pixel locations here. Assume pixels are centered at the integer coordinates.
 * - Handle the image boundaries (coordinates outside [0...w-1] x [0...h-1]) appropriately to avoid artifacts, don't assume
 *   everything outside the boundaries is black or white. One simple method is to just clamp wayward coordinates to the feasible
 *   range.
 * - Make sure each output value is clamped to [0...255] to avoid overflow. Remember to clamp BEFORE type-casting to unsigned
 *   char.
 * - Remember that the image may not have 4 channels. In this case missing channels should be set to 0 in \a sampled_color.
 *
 *
 *   P ------- P ------- P
 *   |         |         |
 *   |  L      |         |
 *   |         |         |
 *   P ------- P ------- P
 *   |         |         |
 *   |         |         |
 *   |         |         |
 *   P ------- P ------- P
 */
void
sampleBilinear(Image const & image, Vec2 const & loc, unsigned char sampled_color[4])
{
  // TODO: REPLACE THIS CODE WITH A CORRECT VERSION

  sampled_color[0] = 0;
  sampled_color[1] = 0;
  sampled_color[2] = 0;
  sampled_color[3] = 0;
}

/**
 * Distorts an image according to the algorithm described in Feature-Based Image Metamorphosis. Linearly interpolates the
 * segments from seg1_start to seg1_end.
 */
Image
distortImage(Image const & image,
             std::vector<LineSegment> const & seg_start,
             std::vector<LineSegment> const & seg_end,
             double t,
             double a, double b, double p)
{
  assert(seg_start.size() == seg_end.size());

  std::cout << "Distorting image..." << std::endl;

  Image result(image.width(), image.height(), image.numChannels());

  // TODO: YOUR CODE HERE

  // Remember to use bilinear interpolation to get the colors of pixels with real-valued row and column coordinates!

  return result;
}

/* Linearly blends corresponding pixels of two images to produce the resulting image. */
Image
blendImages(Image const & img1, Image const & img2, double t)
{
  assert(img1.hasSameDimsAs(img2));

  std::cout << "Blending images..." << std::endl;

  Image result(img1.width(), img1.height(), img1.numChannels());

  // TODO: YOUR CODE HERE

  return result;
}

/* Morph img1 into img2. */
Image
morphImages(Image const & img1,
            Image const & img2,
            std::vector<LineSegment> const & seg1,
            std::vector<LineSegment> const & seg2,
            double t,
            double a, double b, double p)
{
  assert(img1.hasSameDimsAs(img2));

  // TODO: YOUR CODE HERE

  // First distort img1 from 0 to t
  // using seg1 as the initial segments and seg2 as the final ones.
  //
  //    Image distorted1 = ...

  // Then distort img2 from 1 to (1 - t)
  // using seg2 as the initial segments and seg1 as the final ones.
  //
  //    Image distorted2 = ...

  // Now blend the results by linearly interpolating ("lerping")
  //
  //    Image blended = ...

  return Image();  // update this line
}

///////////////////////////////////////////////////////////////////////////////
//
//  Driver functions follow. These should not be modified.
//
///////////////////////////////////////////////////////////////////////////////

/**
 * Read segments defining the map between two images from a text file. Each line of the file consists of a single pair of
 * segments. A segment consists of two 2D points (x, y) defining its start and end. The two segments in a pair identify matching
 * features in the two images.
 */
bool
loadSegments(std::string const & path, std::vector<LineSegment> & seg1, std::vector<LineSegment> & seg2)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    std::cerr << "Could not open correspondence file " << path << std::endl;
    return false;
  }

  seg1.clear();
  seg2.clear();

  double asx, asy, aex, aey, bsx, bsy, bex, bey;
  long num_segs = 0;

  std::string line;
  bool first_line = true;
  while ((first_line || (long)seg1.size() < num_segs) && std::getline(in, line))
  {
    std::istringstream line_in(line);
    if (first_line)
    {
      if (!(line_in >> num_segs))
      {
        std::cerr << "Could not read number of segments";
        return false;
      }

      first_line = false;
    }
    else
    {
      if (!(line_in >> asx >> asy >> aex >> aey >> bsx >> bsy >> bex >> bey))
      {
        std::cerr << "Could not read segment pair " << seg1.size();
        return false;
      }

      seg1.push_back(LineSegment(Vec2(asx, asy), Vec2(aex, aey)));
      seg2.push_back(LineSegment(Vec2(bsx, bsy), Vec2(bex, bey)));
    }
  }

  assert(seg1.size() == seg2.size());

  return (long)seg1.size() == num_segs;
}

bool
morphDriver(std::string const & img1_path, std::string const & img2_path, std::string const & seg_path, double t,
            std::string const & out_path, double a, double b, double p)
{
  // Load images, forcing both to 4-channel RGBA for compatibility
  Image img1, img2;
  if (!img1.load(img1_path, 4) || !img2.load(img2_path, 4))
    return false;

  if (!img1.hasSameDimsAs(img2))
  {
    std::cerr << "Both input images must be the same dimensions" << std::endl;
    return false;
  }

  std::cout << "Loaded two " << img1.width() << 'x' << img1.height() << ' ' << img1.numChannels() << "-channel images"
            << std::endl;

  // Load segments
  std::vector<LineSegment> seg1, seg2;
  if (!loadSegments(seg_path, seg1, seg2))
    return false;

  std::cout << "Read " << seg1.size() << " segments" << std::endl;
  Image morphed = morphImages(img1, img2, seg1, seg2, t, a, b, p);
  if (!morphed.save(out_path))
    return false;

  return true;
}

int
main(int argc, char * argv[])
{
  if (argc != 6 && argc != 9)
  {
    std::cout << "Usage: " << argv[0] << " image1 image2 segments_file time[0..1] output.png [a  b  p]" << std::endl;
    return -1;
  }

  std::string img1_path  =  argv[1];
  std::string img2_path  =  argv[2];
  std::string seg_path   =  argv[3];
  double t               =  std::atof(argv[4]);
  std::string out_path   =  argv[5];

  if (t < 0.0 || t > 1.0)
  {
    std::cout << "Time t out of range: clamping to [0..1]" << std::endl;
    if (t > 1.0)
      t = 1.0;
    else
      t = 0.0;
  }

  double a = 0.5;
  double b = 1;
  double p = 0.2;
  if (argc == 9)
  {
    a = std::atof(argv[6]);
    b = std::atof(argv[7]);
    p = std::atof(argv[8]);
  }

  std::cout << "Morphing " << img1_path << " into " << img2_path << " at time t = " << t << ", generating " << out_path
            << std::endl;
  std::cout << "Using parameters { a : " << a << ", b : " << b << ", p : " << p << " }" << std::endl;

  morphDriver(img1_path, img2_path, seg_path, t, out_path, a, b, p);

  return 0;
}
