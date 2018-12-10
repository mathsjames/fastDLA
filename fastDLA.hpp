#include <vector>
#include <random>
#include <complex>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const int minLength = 4) : root_(nullptr)
  {
    minSideLength = minLength;
    maxDepth_ = 1 + ceil_log2(maxRadius / minSideLength);

    sideLength=new double[maxDepth_];
    double length = minSideLength;
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
  
  struct NodePointers
  {
    Node* NW;
    Node* NE;
    Node* SW;
    Node* SE;
  };

  union Node
  {
    NodePointers pointers;
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
	do
	  {
	    //find direction
	    if (!node->pointers[][])
	      {
		if (depth<scaleDepth && cabs(currPoint)>startDist)
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
		diffuseRecursive(node->pointers[][],depth+1,centre+);
	      }
	  }
	while (particleFree && particleIsPresent());
      }
    else
      {
	do
	  {
	    finishingStep(getNearestPoints(node));
	  }
	while (particleFree && particleIsPresent());
      }
  }

  Node* root_;
  const int maxDepth;
  const double minSideLength;
  double startDist;
  complex<double> currPoint;
  bool particleFree;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif2PI(0.0,2*PI);
}
