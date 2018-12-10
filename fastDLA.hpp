#include <vector>
#include <random>
#include <complex>
#include <cmath>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const int minLength = 4) : root_(nullptr)
  {
    maxDepth_ = 1 + ceil_log2(maxRadius / minLength);

    sideLengths=new double[maxDepth_];
    double length = minLength;
    for (int i = maxDepth_ - 1; i >= 0; i--)
      {
	sideLength[i] = length;
	length*=2;
      }

    currPoint = std::complex<double>(0.0,0.0);
    insertPoint();

    startDist = 2;
  }

  void grow(int n)
  {
    int i;
    for (i = 0; i < n; i++)
      {
	aggregate();
      }
  }

  void writePoints(FILE* fp)
  {
    writeRecursive(root_, fp);
  }

private:

  union Node
  {
    Node pointers[4];
    std::vector<complex<double>> points;
  };

  complex randcirc()
  {
    return exp(std::complex<double>(0,unif2PI(generator)));
  }

  void aggregate()
  {
    currPoint = startdist*randcirc();
    particleFree = true;
    diffuseRecursive(root_,0,std::complex<int>(0,0));
    while (particleFree)
      {
	resetParticle();
	diffuseRecursive(root_,0,std::complex<int>(0,0));
      }
    insertPoint();
  }

  void insertPoint()
  {
    startDist=std::max(startDist,cabs(currPoint)+2);
    
    if (startDist>maxRadius)
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw Exception();
      }
    
    if (startDist>sideLength[scaleDepth])
      scaleDepth--;
    
    insertRecursive(root_, 0, std::complex<int>(0,0));
  }

  void insertRecursive(Node* node, int depth, complex<int> centre)
  {
  }

  struct Nearest
  {
    complex<double> nearest;
    double nextDist;
  };

  void diffuseRecursive(Node* node, int depth, complex<int> centre)
  {
    if (depth<maxDepth)
      {
	int dirx, diry, dir;
	do
	  {
	    dirx = (real(currPoint)>real(centre));
	    diry = (imag(currPoint)>imag(centre));
	    dir = dirx+2*diry;
	    if (!node->pointers[dir])
	      {
		if (depth<=scaleDepth && cabs(currPoint)>startDist)
		  {
		    resetParticle();
		  }
		else
		  {
		    step(sideLength[depth]);
		  }
	      }
	    else
	      {
		int nextCentre = centre+(sideLength[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
		diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
	      }
	  }
	while (particleFree && particleIsPresent(depth, centre));
      }
    else
      {
	do
	  {
	    finishingStep(getNearestPoints(node, centre));
	  }
	while (particleFree && particleIsPresent(depth, centre));
      }
  }

  void step(double size)
  {
    currPoint+=size*randCirc();
  }

  void resetParticle()
  {
  }

  void finishingStep(Nearest nearestInfo)
  {
  }

  void getNearestPoints(Node* node, complex<int> centre)
  {
  }

  bool particleIsPresent(int depth, complex<int> centre)
  {
    complex<double> disp = real(currPoint);
    double icp=imag(currPoint);
    return ( abs(rcp-real(centre))<sideLengths[depth] &&
	     abs(icp-imag(centre))<sideLengths[depth] );
  }

  Node* root_;
  const int maxDepth;
  int scaleDepth;
  int sideLengths[];
  double startDist;
  complex<double> currPoint;
  bool particleFree;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif2PI(0.0,2*PI);
}
