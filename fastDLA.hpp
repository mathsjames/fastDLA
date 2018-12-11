#include <vector>
#include <random>
#include <complex>
#include <cmath>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const int halfMinLength = 2) : root_(nullptr)
  {
    maxDepth_ = 1 + ceil(log2(maxRadius / minLength));

    sideLengths=new int[maxDepth_+1];
    int length = halfMinLength;
    for (int i = maxDepth_; i >= 0; i--)
      {
	sideLength[i] = length;
	length*=2;
      }

    currPoint = std::complex<double>(0.0,0.0);
    insertPoint();

    startDist = 2;
  }

  ~ClusterTree()
  {
    clearRecursive();
    delete sideLengths;
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
    Node* pointers[4];
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
    double distToNearest;
    double maxSafeDist;
  };

  struct gridEntry
  {
    complex<int> centre;
    double minDist;
  }

  void diffuseRecursive(Node* node, int depth, complex<int> centre)
  {
    while (particleFree)
      {
	if (!gettingNearest)
	  {
	    if (!particleIsPresent(depth, centre))
	      break;

	    if (depth!=maxDepth)
	      {
		int dirx, diry, dir;
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
			step(sideLength[depth]-2);
		      }
		  }
		else
		  {
		    int nextCentre = centre+(sideLength[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
		    diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
		  }
	      }
	    else
	      {
		gettingNearest = true;
		nearestInfo.maxSafeDist2=3*sideLengths[depth]-std::max(abs(real(currPoint-centre)),abs(imag(currPoint-centre)));
		nearestInfo.maxSafeDist2=nearestInfo.maxSafeDist2*nearestInfo.maxSafeDist2;
		nearestInfo.distToNearest2=nearestInfo.maxSafeDist2;
		initialiseGrid(centre);
	      }
	  }
	else
	  {
	    if (grid.size()==0)
	      {
		if (nearestInfo.dist2ToNearest>=nearestInfo.maxSafeDist2)
		  {
		    step(nearestInfo.maxSafeDist);
		  }
		else
		  {
		    finishingStep(nearestInfo);
		  }
		gettingNearest = false;
	      }
	    else
	      {
		if (depth==maxDepth_)
		  {
		    checkForCloser(node,centre);
		    break;
		  }
		else
		  {
		    int dirx, diry;
		    if (gridBeneath(centre,depth,dirx,diry))
		      {
			int dir = dirx+2*diry;
			if (node->pointers[dir])
			  {
			    int nextCentre = centre+(sideLength[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
			    diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
			  }
			else
			  {
			    markEmpty(centre,depth);
			  }
		      }
		    else
		      {
			break;
		      }
		  }
	      }
	  }
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

  void initializeGrid(complex<int> centre)
  {
    int i, j;
    for (i = -1; i < 2; i++)
      {
	for (j = -1; j < 2; j++)
	  {
	    gridEntry entry;
	    entry.centre = centre+2*sideLength[maxDepth_]*std::complex<int>(i,j);
	    entry.dist2 = squaredDistToSquare(centre);
	    grid.push_back(entry);
	  }
      }
  }

  double squaredDistToSquare(complex<int> centre)
  {
    double dist2 = 0;
    double xOffset = abs(real(centre)-real(currPoint));
    double yOffset = abs(imag(centre)-imag(currPoint));
    if (x>sideLength[maxDepth_]) dist2+=(x-sideLength[maxDepth_])*(x-sideLength[maxDepth_]);
    if (y>sideLength[maxDepth_]) dist2+=(y-sideLength[maxDepth_])*(y-sideLength[maxDepth_]);
    return dist2;
  }

  void checkForCloser(Node* node, complex<int> centre)
  {
    for (int i = 0; i < node->points.size(); i++)
      {
	complex<double> newPoint = node->points[i];
	double dist2 = squaredDist(newPoint);
	if (dist2<nearestInfo.maxSafeDist2)
	  {
	    if (dist2<nearestInfo.dist2ToNearest)
	      {
	      nearestInfo.maxSafeDist2 = nearestInfo.dist2ToNearest;
	      nearestInfo.dist2ToNearest = dist2;
	      nearestInfo.nearest = newPoint;
	      }
	    else
	      {
		nearestInfo.maxSafeDist2 = dist2;
	      }
	  }
      }
    updateGrid(centre);
  }

  void updateGrid(complex<int> centre)
  {
    int i;
    for (i = 0; i < grid.size(); i++)
      {
	if (grid[i].centre=centre || grid[i].dist2>nearestInfo.maxSafeDist2)
	  {
	    //remove from vector
	  }
      }
  }

  void markEmpty(complex<int> centre,int depth)
  {
    int i;
    for (i = 0; i < grid.size(); i++)
      {
	complex<int> diff = grid[i].centre-centre;
	if ( abs(real(diff))<sideLengths[depth] && abs(imag(diff))<sideLengths[depth])
	  {
	    //remove from vector
	  }
      }
  }

  bool gridBeneath(complex<int> centre,int depth, int& dirx, int& diry)
  {
  }

  bool particleIsPresent(int depth, complex<int> centre)
  {
    complex<double> disp = real(currPoint);
    return ( abs(real(disp))<sideLengths[depth] &&
	     abs(imag(disp))<sideLengths[depth] );
  }

  Node* root_;
  const int maxDepth;
  int scaleDepth;
  int sideLengths[];
  double startDist;
  complex<double> currPoint;
  bool particleFree;
  Nearest nearestInfo;
  std::vector<gridEntry> grid;
  bool gridEmpty;
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif2PI(0.0,2*PI);
}
