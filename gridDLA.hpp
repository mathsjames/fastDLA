#include <vector>
#include <random>
#include <chrono>
#include <complex>
#include <exception>
#include <cmath>

#ifndef SECOND_H
#define SECOND_H
#include "baseDLA.hpp"
#endif

class ClusterGrid : public baseDLA
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
    
    pointsGridSize = floor(2*maxRadius/minPointGridMesh);
    pointsGridMesh = 2*maxRadius/pointsGridSize;
    pointsGrid = new std::vector<int>[pointsGridSize*pointsGridSize];
    std::fill(pointsGrid, pointsGrid+pointsGridSize*pointsGridSize, std::vector<int>(0));

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

  ~ClusterGrid()
  {
    delete layerSizes;
    delete layerMeshes;
    for (int i=0; i<layerCount; i++)
      {
	delete (layers[i]);
      }
    delete layers;
    delete pointsGrid;
    delete points;
  }

  void writePoints(FILE* fp)
  {
    fwrite(points,sizeof(std::complex<double>),pointsAdded,fp);
  }

private:
  void grow(int n)
  {
    int i;
    for (i = 0; i < n; i++)
      {
	aggregate();
      }
  }

  void aggregate()
  {
    setStartPoint();
    particleFree = true;
    while (particleFree)
      {
	findAndMakeStep();
      }
    addParticle();
    markParticle();
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
  };
    
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
		nearestInfo.nearest=points[vec[k]];
		nearestInfo.maxSafeDist2 = nearestInfo.distToNearest2;
		nearestInfo.distToNearest2 = dist2;
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
  }

  void step(int layer)
  {
    currPoint += ((double) layerSizes[layer])*randCirc();
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
    pointsGrid[indices.index1*pointsGridSize+indices.index2].push_back(pointsAdded);
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
    for (int layer = layerCount - 1; layer >= 0; layer--)
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
    
    if (startDist>maxRadius_)
      {
	fprintf(stderr,"Cluster Exceeded maxRadius\n");
	throw 373;
      }
  }

  double maxRadius_;
  double startDist;
  
  int layerCount;
  int* layerSizes;
  double* layerMeshes;
  char** layers;
  
  int pointsGridSize;
  double pointsGridMesh;
  std::vector<int>* pointsGrid;
  
  std::complex<double>* points;
  int pointsAdded;

  std::complex<double> currPoint;
  bool particleFree;

  std::mt19937 generator;
  std::uniform_real_distribution<double> unif2PI;
  std::cauchy_distribution<double> cauchy;

  const std::complex<double> cx_1 = std::complex<double>(1.0,0.0);
  const double PI = 3.1415926535;
  const double two_PI = 6.2831853072;
};
