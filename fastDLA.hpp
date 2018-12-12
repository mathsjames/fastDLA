#include <vector>
#include <random>
#include <complex>
#include <cmath>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const int halfMinLength = 2) : root_(nullptr)
  {
    maxDepth_ = ceil(log2(maxRadius / halfMinLength));

    sideLengths = malloc(sizeof(int)*(maxDepth_+1));
    int length = halfMinLength;
    for (int i = maxDepth_; i >= 0; i--)
      {
	sideLengths[i] = length;
	length*=2;
      }

    root_ = new Node();

    currPoint = std::complex<double>(0.0,0.0);
    Node* node = root_->pointers[3] = new Node();
    int depth=1;
    while (depth<maxDepth_)
      {
	node = node->pointers[0] = new Node();
	depth++;
      }
    node->points.push_back(currPoint);
    startDist = 2;
    markPoint();
  }

  ~ClusterTree()
  {
    clearRecursive(root_,depth);
    free(sideLengths);
    delete root_;
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
    std::vector<std::complex<double>> points;
  };

  void clearRecursive(Node* node, int depth)
  {
    if (depth!=maxDepth)
      {
	int i;
	for (i=0;i<4;i++)
	  {
	    if (node->pointers[i])
	      {
		clearRecursive(node->pointers[i],depth+1);
		delete node->pointers[i];
	      }
	  }
      }
  }

  std::complex<double> randCirc()
  {
    return exp(std::complex<double>(0,unif2PI(generator)));
  }

  void aggregate()
  {
    currPoint = startDist*randCirc();
    particleFree = true;
    needToMark = true;
    diffuseRecursive(root_,0,std::complex<int>(0,0));
    while (particleFree)
      {
	resetParticle();
	diffuseRecursive(root_,0,std::complex<int>(0,0));
      }
    updateStartDist();
    if (needToMark)
      {
	markPoint();
      }
  }

  void updateStartDist()
  {
    startDist=std::max(startDist,cabs(currPoint)+2);
    
    if (startDist>sideLengths[0])
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw Exception();
      }
    
    if (startDist>sideLengths[scaleDepth])
      scaleDepth--;
  }

  void markPoint()
  {
    markRecursive(root_, 0, std::complex<int>(0,0));
  }

  void markRecursive(Node* node, int depth, std::complex<int> centre)
  {
    if (depth!=maxDepth_)
      {
	int dirx, diry, dir;
	for (dirx = 0; dirx < 2; dirx++)
	  {
	    for (diry = 0; diry < 2; diry++)
	      {
		std::complex<int> nextCentre = centre+sideLengths[depth+1]*std::complex<int>(2*dirx-1,2*diry-1);
		if (isWithinLInf(currPoint,nextCentre,3*sideLengths[depth+1]))
		  {
		    dir = dirx+2*diry;
		    if (!node->pointers[dir])
		      {
			node->pointers[dir] = new Node();
		      }
		    markRecursive(node->pointers[dir],depth+1,nextCentre);
		  }
	      }
	  }
      }
  }

  bool isWithinLInf(std::complex<double> point1,std::complex<double> point2,double dist)
  {
    std::complex<double> disp = point1-point2;
    return ( abs(real(disp))<=dist &&
	     abs(imag(disp))<=dist );
  }

  struct Nearest
  {
    std::complex<double> nearest;
    double distToNearest2;
    double maxSafeDist2;
  };

  struct gridEntry
  {
    std::complex<int> centre;
    double minDist2;
  };

  void diffuseRecursive(Node* node, int depth, std::complex<int> centre)
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
			step(sideLengths[depth]-2);
		      }
		  }
		else
		  {
		    int nextCentre = centre+(sideLengths[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
		    diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
		  }
	      }
	    else
	      {
		gettingNearest = true;
		nearestInfo.maxSafeDist2=3*sideLengths[depth]-std::max(abs(real(currPoint-centre)),abs(imag(currPoint-centre)));
		nearestInfo.maxSafeDist2=nearestInfo.maxSafeDist2*nearestInfo.maxSafeDist2;
		nearestInfo.distToNearest2=nearestInfo.maxSafeDist2;
		initializeGrid(centre);
	      }
	  }
	else
	  {
	    if (grid.size()==0)
	      {
		if (nearestInfo.distToNearest2>=nearestInfo.maxSafeDist2)
		  {
		    step(sqrt(nearestInfo.maxSafeDist2));
		  }
		else
		  {
		    finishingStep(nearestInfo);
		    if (!particleFree)
		      {
			if (!node->points.empty())
			  {
			    needToMark = false;
			  }
			node->points.push_back(currPoint);
		      }
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

  void initializeGrid(std::complex<int> centre)
  {
    int i, j;
    for (i = -1; i < 2; i++)
      {
	for (j = -1; j < 2; j++)
	  {
	    gridEntry entry;
	    entry.centre = centre+2*sideLengths[maxDepth_]*std::complex<int>(i,j);
	    entry.minDist2 = squaredDistToSquare(centre);
	    grid.push_back(entry);
	  }
      }
  }

  double squaredDistToSquare(std::complex<int> centre)
  {
    double dist2 = 0;
    double x = abs(real(centre)-real(currPoint));
    double y = abs(imag(centre)-imag(currPoint));
    if (x>sideLength[maxDepth_]) dist2+=(x-sideLengths[maxDepth_])*(x-sideLengths[maxDepth_]);
    if (y>sideLength[maxDepth_]) dist2+=(y-sideLengths[maxDepth_])*(y-sideLengths[maxDepth_]);
    return dist2;
  }

  void checkForCloser(Node* node, std::complex<int> centre)
  {
    for (int i = 0; i < node->points.size(); i++)
      {
	std::complex<double> newPoint = node->points[i];
	double dist2 = squaredDist(currPoint,newPoint);
	if (dist2<nearestInfo.maxSafeDist2)
	  {
	    if (dist2<nearestInfo.distToNearest2)
	      {
		nearestInfo.maxSafeDist2 = nearestInfo.distToNearest2;
		nearestInfo.distToNearest2 = dist2;
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

  double squaredDist(std::complex<double> point1, std::complex<double> point2)
  {
    std::complex<double> disp = point1-point2;
    return real(disp)*real(disp)+imag(disp)*imag(disp);
  }

  void updateGrid(std::complex<int> centre)
  {
    grid.remove_if( [&](gridEntry entry)
		    {
		      return entry.centre==centre || entry.minDist2>nearestInfo.maxSafeDist2;
		    });
  }

  void markEmpty(std::complex<int> centre,int depth)
  {
    grid.remove_if( [&](gridEntry entry)
		    {
		      std::complex<int> diff = entry.centre-centre;
		      return (abs(real(diff))<sideLengths[depth] && abs(imag(diff))<sideLengths[depth]);
		    });
}

  bool gridBeneath(std::complex<int> centre,int depth, int& dirx, int& diry)
  {
    bool retval = false;

    std::list<gridEntry>::iterator it;
    gridEntry entry;

    for (it=grid.begin(); it != grid.end(); ++it)
      {
	entry=*it;
	std::complex<int> diff = entry.centre-centre;
	if ( abs(real(diff))<sideLengths[depth] && abs(imag(diff))<sideLengths[depth])
	  {
	    dirx = sign(real(diff));
	    diry = sign(imag(diff));
	    retval = true;
	    break;
	  }
      }
    return retval;
  }

int sign(int val)
{
  if (val>0) return 1;
  if (val<0) return -1;
  return 0;
}

bool particleIsPresent(int depth, std::complex<int> centre)
{
  return isWithinLInf(currPoint,centre,sideLengths[depth]);
}

Node* root_;
const int maxDepth_;
int scaleDepth;
int sideLengths[];
double startDist;
std::complex<double> currPoint;
bool particleFree;
bool gettingNearest;
Nearest nearestInfo;
std::list<gridEntry> grid;
bool gridEmpty;
bool needToMark;
std::default_random_engine generator;
std::uniform_real_distribution<double> unif2PI(0.0,2*PI);
};
