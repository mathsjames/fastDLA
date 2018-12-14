#include "fastDLA.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  int n, seed;
  char filename[1000];

  seed=time(NULL);
  filename[0]='\0';

  switch(argc-1)
    {
    case 3:
      sscanf(argv[3],"%d",&seed);
    case 2:
      sscanf(argv[2],"%s",&filename);
    case 1:
      sscanf(argv[1],"%d",&n);
      break;
    default:
      fprintf(stderr,"First argument (Cluster Size) is mandatory.\nSecond argument (Filename to save to) and third argument (seed) are optional.\n");
      return 1;
    }

  double maxRadius = 5.0+2*pow(n,0.75); //for asymptically big enough maxRadius need exponent to be bigger than 1/d where d is dimension of DLA clusters
  ClusterTree cluster(maxRadius);

  std::cout << "Starting to grow cluster" << std::endl;
  cluster.grow(n-1);//one particle beside the root was already added at construction
  std::cout << "Finished growing cluster" << std::endl;

  if (filename[0])
    {
      FILE* fp=fopen(filename, "w");
      if (!fp)
	{
	  fprintf(stderr, "Failed to open file to save\n");
	  return 1;
	}
      else
	{
	  cluster.writePoints(fp);
	  fclose(fp);
	  printf("Written cluster to file\n");
	}
    }
  else
    {
      printf("Not saving points as no filename given\n");
    }

  return 0;
}
