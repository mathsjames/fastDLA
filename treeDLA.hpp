#include <vector>
#include <list>
#include <random>
#include <exception>
#include <complex>
#include <cmath>
#include <chrono>
#include <iostream>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const int halfMinLength = 2, unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count()) : root_(nullptr)
  {
    maxDepth_ = ceil(log2(maxRadius / halfMinLength));
    //std::cout << "maxDepth" << maxDepth_ << std::endl;
    //std::cout << "sizeof Node" << sizeof(Node) << std::endl;
    sideLengths =(int*) malloc(sizeof(int)*(maxDepth_+1));
    int length = halfMinLength;
    for (int i = maxDepth_; i >= 0; i--)
      {
	sideLengths[i] = length;
	length*=2;
      }

    generator.seed(seed);
    unif2PI = std::uniform_real_distribution<double>(0.0,two_PI);
    cauchy = std::cauchy_distribution<double>(0.0,1.0);

    root_ = new Node();

    std::complex<double> zeroPoint = std::complex<double>(0.0,0.0);
    std::complex<double> firstPoint = 2.0*randCirc();
    int dir = (std::real(firstPoint)>0)+2*(std::imag(firstPoint)>0);
    Node* node = root_->pointers[dir] = new Node();
    int depth=1;
    while (depth<maxDepth_)
      {
	node = node->pointers[3-dir] = new Node();
	depth++;
      }
    node->points.push_back(zeroPoint);
    node->points.push_back(firstPoint);
    //std::cout << "First Point " << firstPoint << std::endl;
    currPoint = firstPoint;
    markPoint();

    startDist = 4;
    scaleDepth = maxDepth_;

    gettingNearest = false;
    addingPoint = false;
  }

  ~ClusterTree()
  {
    clearRecursive(root_, 0);
    free(sideLengths);
    delete root_;
  }

  void grow(int n)
  {
    int i;
    for (i = 0; i < n; i++)
      {
	//std::cout << "aggregate " << i << std::endl;
	aggregate();
      }
  }

  void writePoints(FILE* fp)
  {
    writeRecursive(root_, 0, fp);
  }

private:

  struct Node
  {
    Node* pointers[4];
    std::vector<std::complex<double> > points;
  };

  void clearRecursive(Node* node, int depth)
  {
    if (depth!=maxDepth_)
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

  void writeRecursive(Node* node, int depth, FILE* fp)
  {
    if (depth!=maxDepth_)
      {
	int i;
	for (i=0;i<4;i++)
	  {
	    if (node->pointers[i])
	      {
		//std::cout << "going down to depth " << depth << " in direction " << i << std::endl;
		writeRecursive(node->pointers[i], depth+1, fp);
		//std::cout << "going up" << std::endl;
	      }
	  }
      }
    else
      {
	//std::cout << "writing " << node->points.size() << " from node " << node << std::endl;
	fwrite(&(node->points[0]),sizeof(std::complex<double>),node->points.size(),fp);
      }
  }

  std::complex<double> randCirc()
  {
    return std::exp(std::complex<double>(0,unif2PI(generator)));
  }

  void aggregate()
  {
    currPoint = startDist*randCirc();
    //std::cout << "Starting diffusing at " << currPoint << std::endl;
    particleFree = true;
    needToMark = true;
    diffuseRecursive(root_,0,std::complex<int>(0,0));
    if (particleFree)
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw 373;
      }
    //std::cout << "Done diffusing at " << std::real(currPoint) << "+i" << std::imag(currPoint) << std::endl;
    updateStartDist();
    //std::cout << "Start radius is now " << startDist << std::endl;
    if (needToMark)
      {
	//std::cout << "Marking" << std::endl;
	markPoint();
      }
  }

  void updateStartDist()
  {
    startDist=std::max(startDist,abs(currPoint)+2);
    
    if (startDist>sideLengths[0])
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw 373;
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
		if (isWithinLInf(currPoint,cast(nextCentre),3*sideLengths[depth+1]))
		  {
		    dir = dirx+2*diry;
		    if (!node->pointers[dir])
		      {
			//std::cout << "creating node with centre " << nextCentre << std::endl;
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
    return ( abs(std::real(disp))<=dist &&
	     abs(std::imag(disp))<=dist );
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
	//std::cout << "currPoint " << currPoint << std::endl;
	//std::cout << "depth:" << depth << " node:" << (long) node << " centre:" << real(centre) << "+i" << imag(centre) << std::endl;
	if (!gettingNearest)
	  {
	    //std::cout << "looking for particle" << std::endl;
	    if (!particleIsPresent(depth, centre))
	      break;
	    //std::cout << "particle present" << std::endl;
	    //std::cout << "depth:" << depth << "  maxDepth:" << maxDepth_ << std::endl;
	    if (depth!=maxDepth_)
	      {
		int dirx, diry, dir;
		dirx = (std::real(currPoint)>std::real(centre));
		diry = (std::imag(currPoint)>std::imag(centre));
		dir = dirx+2*diry;
		if (!node->pointers[dir])
		  {
		    if (depth<=scaleDepth && abs(currPoint)>startDist)
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
		    std::complex<int> nextCentre = centre+(sideLengths[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
		    diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
		  }
	      }
	    else
	      {
		gettingNearest = true;
		nearestInfo.maxSafeDist2=3*sideLengths[depth]-std::max(abs(std::real(currPoint-cast(centre))),abs(std::imag(currPoint-cast(centre))));
		nearestInfo.maxSafeDist2=nearestInfo.maxSafeDist2*nearestInfo.maxSafeDist2;
		nearestInfo.distToNearest2=nearestInfo.maxSafeDist2;
		initializeGrid(centre);
	      }
	  }
	else
	  {
	    //std::cout << "looking for nearby particles" << std::endl;
	    //std::cout << "grid size:" << grid.size() << std::endl;
	    if (grid.size()==0)
	      {
		if (nearestInfo.distToNearest2>=nearestInfo.maxSafeDist2)
		  {
		    //std::cout << "distToNearest2 " << nearestInfo.distToNearest2 << std::endl;
		    //std::cout << "maxSafeDist2 " << nearestInfo.maxSafeDist2 << std::endl;
		    step(sqrt(nearestInfo.maxSafeDist2)-2);
		  }
		else
		  {
		    finishingStep();
		    if (!particleFree)
		      {
			addingPoint = true;
		      }
		  }
		gettingNearest = false;
	      }
	    else
	      {
		//std::cout << "depth:" << depth << "  maxDepth:" << maxDepth_ << std::endl;
		if (depth==maxDepth_)
		  {
		    checkForCloser(node,centre);
		    break;
		  }
		else
		  {
		    int dirx, diry; // passed by reference to be set bu gridBeneath
		    if (gridBeneath(centre,depth,dirx,diry))
		      {
			int dir = dirx+2*diry;
			//std::cout << "dir (x,y,):" << dirx << diry <<dir << std::endl;
			if (node->pointers[dir])
			  {
			    std::complex<int> nextCentre = centre+(sideLengths[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
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
    if (addingPoint && particleIsPresent(depth, centre))
      {
	if (depth==maxDepth_)
	  {
	    if (!node->points.empty())
	      {
		needToMark = false;
	      }
	    //std::cout << "adding particle at " << currPoint << " to node " << node << " with depth " << depth << std::endl;
	    node->points.push_back(currPoint);
	    addingPoint=false;
	  }
	else
	  {
	    int dirx, diry, dir;
	    dirx = (std::real(currPoint)>std::real(centre));
	    diry = (std::imag(currPoint)>std::imag(centre));
	    dir = dirx+2*diry;
	    std::complex<int> nextCentre = centre+(sideLengths[depth]/2)*std::complex<int>(2*dirx-1,2*diry-1);
	    diffuseRecursive(node->pointers[dir],depth+1,nextCentre);
	  }
      }
  }

  void step(double size)
  {
    //std::cout << "stepping from " << currPoint << " to ";
    currPoint+=size*randCirc();
    //std::cout << currPoint << std::endl;
  }

  void resetParticle()
  {
    //std::cout << "resetting particle" << std::endl;
    double ratioOut = abs(currPoint)/startDist;
    std::complex<double> z = harmonicToCircle(ratioOut);
    currPoint *= z/ratioOut;
  }
  
  std::complex<double> harmonicToCircle(double absPos)
  {
    std::complex<double> z = randCirc();
    z = (z-cx_1)/(z+cx_1);
    z *= (absPos-1)/(absPos+1);
    z = -(z+cx_1)/(z-cx_1);
    return z;
  }

  void finishingStep()
  {
    //std::cout << "fancy stepping from " << currPoint << std::endl;
    //std::cout << "Nearest:" << nearestInfo.nearest << std::endl;
    //std::cout << "squared distance:" << nearestInfo.distToNearest2 << std::endl;
    //std::cout << "squared next distance:" << nearestInfo.maxSafeDist2 << std::endl;
    double d1 = sqrt(nearestInfo.distToNearest2)/2;
    double d2 = sqrt(nearestInfo.maxSafeDist2)/2;
    double theta=acos((1+(d2-1)*(d2-1)-d1*d1)/(2*(d2-1)));
    double r2=((d2-1)*(d2-1)+d1*d1-1)/(2*d1);
    double r1=sqrt(1-(d1-r2)*(d1-r2));
    std::complex<double> alpha=std::complex<double>(r1,-r2);
    std::complex<double> beta=std::complex<double>(-r1,-r2);
    std::complex<double> D=(d2-1-beta)/(d2-1-alpha);
    std::complex<double> y1=std::pow(alpha*D/beta,PI/theta);
    double y2=real(y1)+imag(y1)*cauchy(generator);
    std::complex<double> y3=pow(std::complex<double>(y2,0.0),theta/PI);
    std::complex<double> y4=(-beta*y3+D*alpha)/(-y3+D);
    currPoint += (std::complex<double>(0.0,1.0)*y4*(nearestInfo.nearest-currPoint)/d1);
    particleFree = (y2>=0);
    //std::cout << "stepped to " << currPoint << std::endl;
    //std::cout << "finished?:" << !particleFree << std::endl;
    //std::cout << "Report:" << d1 << " " << d2 << " " << theta << " " << r2 << " " << r1 << " " << alpha << " " << beta << " " << D << " " << y1 << " " << y2 << " " << y3 << " " << y4 << std::endl;
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
    double x = abs(std::real(centre)-std::real(currPoint));
    double y = abs(std::imag(centre)-std::imag(currPoint));
    if (x>sideLengths[maxDepth_]) dist2+=(x-sideLengths[maxDepth_])*(x-sideLengths[maxDepth_]);
    if (y>sideLengths[maxDepth_]) dist2+=(y-sideLengths[maxDepth_])*(y-sideLengths[maxDepth_]);
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
    return std::real(disp)*std::real(disp)+std::imag(disp)*std::imag(disp);
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
		      return (abs(std::real(diff))<sideLengths[depth] && abs(std::imag(diff))<sideLengths[depth]);
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
	if ( abs(std::real(diff))<sideLengths[depth] && abs(std::imag(diff))<sideLengths[depth])
	  {
	    dirx = (std::real(diff)>=0);
	    diry = (std::imag(diff)>=0);
	    retval = true;
	    break;
	  }
      }
    return retval;
  }

  bool particleIsPresent(int depth, std::complex<int> centre)
  {
    return isWithinLInf(currPoint,cast(centre),sideLengths[depth]);
  }

  std::complex<double> cast(std::complex<int> z)
  {
    return std::complex<double>(std::real(z),std::imag(z));
  }

  Node* root_;
  int maxDepth_;
  int scaleDepth;
  int *sideLengths;
  double startDist;
  std::complex<double> currPoint;
  bool particleFree;
  bool gettingNearest;
  bool addingPoint;
  Nearest nearestInfo;
  std::list<gridEntry> grid;
  bool gridEmpty;
  bool needToMark;
  std::mt19937 generator;
  std::uniform_real_distribution<double> unif2PI;
  std::cauchy_distribution<double> cauchy;

  const std::complex<double> cx_1 = std::complex<double>(1.0,0.0);
  const double PI = 3.1415926535;
  const double two_PI = 6.2831853072;
};
