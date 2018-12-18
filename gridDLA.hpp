#include <vector>
#include <random>
#include <chrono>
#include <complex>
#include <exception>
#include <cmath>

class ClusterGrid
{
public:
  ClusterGrid(const int numberOfParticles, const double maxRadius, const double minLength = 1, const double minPointGridMesh = 4, unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count())
  {
    maxRadius_ = maxRadius;
    
    layerCount = 1 + ceil(log2(maxRadius/minLength));
    layerSizes = new int[layerCount];
    layerMeshes = new double[layerCount];
    layers = new char*[layerCount];
    layerSizes[0]=2;
    for (int i = 0; i < layerCount; i++)
      {
	if (i>0)
	  {
	    layerSizes[i]=layerSizes[i-1]*2;
	  }
	layerMeshes[i] = 2*maxRadius/layerSizes[i];
	layers[i] = new char[layerSizes[i]*layerSizes[i]];
	std::fill(layers[i],layers[i]+layerSizes[i]*layerSizes[i],0);
      }
    
    pointGridSize = floor(2*maxRadius/minPointGridMesh);
    pointGridMesh = 2*maxRadius/pointGridSize;
    pointGrid = new std::vector<int>[pointGridSize*pointGridSize];
    std::fill(pointGrid, pointGrid+pointGridSize*pointGridSize, std::vector<int>(0));

    points = new std::complex<double>[numberOfParticles+1];
    pointsAdded = 0;

    currPoint = std::complex<double>(0.0,0.0);
    addParticle();
    markParticle();
    currPoint = 2.0*randCirc();
    addParticle();
    markParticle();

    startDist = 4;
    
    generator.seed(seed);
    unif2PI = std::uniform_real_distribution<double>(0.0,two_PI);
    cauchy = std::cauchy_distribution<double>(0.0,1.0);

    grow(numberOfParticles-1);
  }

  void grow(int n)
  {
    int i;
    for (i = 0; i < n; i++)
      {
	aggregate();
      }
  }

private:

  void aggregate()
  {
    setStartPoint();
    particleFree = true;
    while (particleFree)
      {
	findAndMakeStep();
      }
    addParticle(currPoint);
    markParticle(currPoint);
    updateStartDist();
  }



  void setStartPoint()
  {
    currPoint = startDist*randCirc();
  }

  std::complex<double> randCirc()
  {
    return std::exp(std::complex<double>(0,unif2PI(generator)));
  }

  struct NearestInfo
  {
    std::complex<double> nearest;
    double distToNearest2;
    double maxSafeDist2;
  };
  
  void findAndMakeStep()
  {
    int best = findBestLayer();
    if (best == layerCount)
      {
	NearestInfo nearestInfo = findNearest();
	finishingStep(nearestInfo);
      }
    else
      {
	step(best);
      }
  }

  int findBestLayer()
  {
    int lowerBound = 0;
    int upperBound = layerCount;
    int midpoint;
    while (lowerBound != upperBound)
      {
	midpoint = (lowerBound + upperBound)/2;
	if (isMarkedAtLayer(midpoint))
	  {
	    lowerBound = midpoint + 1;
	  }
	else
	  {
	    upperBound = midpoint;
	  }
      }
    return lowerBound;
  }

  struct Indices
  {
    int index1;
    int index2;
  }
    
  bool isMarkedAtLayer(int layer)
  {
    Indices indices = getIndices(layer);
    return (bool) layers[layer][indices.index1*layerSizes[layer]+indices.index2];
  }

  NearestInfo findNearest()
  {
    NearestInfo nearestInfo;
    nearestInfo.maxSafeDist2 = pointsGridMesh*pointsGridMesh;
    Indices indices = getIndices(pointsGridMesh);
    int i, j;
    for (i=std::max(indices.index1-1,0);i<std::min(indices.index1+2,pointsGridSize);i++)
      {
	for (j=std::max(indices.index2-1,0);j<std::min(indices.index2+2,pointsGridSize);j++)
	  {
	    checkForCloser(i,j,nearestInfo);
	  }
      }
  }

  void checkForCloser(int i, int j, NearestInfo& nearestInfo)
  {
    std::vector<int> vec = pointsGrid[i*pointsGridSize+j];
    double dist2;
    int k;
    for (k=0;k<vec.size();k++)
      {
	dist2 = abs(currPoint-points[vec[k]]);
	if (dist2<nearestInfo.maxSafeDist2)
	  {
	    if (dist2<nearestInfo.distToNearest2)
	      {
		nearest=points[vec[k]];
		nearestInfo.maxSafeDist2 = nearestInfo.distToNearest2;
		nearestInfo.distToNearest = dist2;
	      }
	    else
	      {
		nearestInfo.maxSafeDist2 = dist2;
	      }
	  }
      }
  }

  
  void finishingStep(NearestInfo nearestInfo)
  {
    
  }

  void step(int layer)
  {
    currPoint += layerSizes[layer]*randCirc();
  }

  void addParticle()
  {
    points[pointsAdded] = currPoint;
    addToGrid();
    pointsAdded++;
  }

  void addToGrid()
  {
    Indices indices = getIndices(pointsGridMesh);
    pointsGrid[indices.index1*pointGridSize+indices.index2].push_back(pointsAdded);
  }

  Indices getIndices(double mesh)
  {
    Indices indices;
    indices.index1 = floor((maxRadius_-std::imag(currPoint))/mesh);
    indices.index2 = floor((maxRadius_+std::real(currPoint))/mesh);
    return indices;
  }

  void markParticle()
  {
    Indices indices;
    for (layer = layerCount - 1; layer >= 0; layer--)
      {
	indices = getIndices(layerMeshes[layer]);
	setMarks(layer,indices);
      }
  }

  void setMarks(int layer, Indices indices)
  {
    int i, j;
    for (i=std::max(indices.index1-1,0);i<std::min(indices.index1+2,layerSizes[layer]);i++)
      {
	for (j=std::max(indices.index2-1,0);j<std::min(indices.index2+2,layerSizes[layer]);j++)
	  {
	    layers[layer][i*layerSizes[layer]+j]=1;
	  }
      }
  }

  void updateStartDist()
  {
    startDist = std::max(startDist, abs(currPoint) + 2);
    
    if (startDist>maxRadius)
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw 373;
      }
  }

  double maxRadius_;
  
  int layerCount;
  int* layerSizes;
  double* layerMeshes;
  char** layers;
  
  int pointGridSize;
  double pointGridMesh;
  std::vector<int>* pointGrid;
  
  std::complex<double>* points;
  int pointsAdded;

  std::complex<double> currPoint;
  bool particleFree

  std::mt19937 generator;
  std::uniform_real_distribution<double> unif2PI;
  std::cauchy_distribution<double> cauchy;
}
