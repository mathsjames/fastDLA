#include "fastDLA.hpp"

int main(int argc, int argv)
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
    }

  double maxRadius = 5.0+pow(n,0.75);
  ClusterTree cluster(maxRadius);

  printf("Starting to grow cluster\n");
  cluster.grow(n);
  printf("Finished growing cluster\n");

  if (filename[0])
    {
      FILE* fp=fopen(filename, 'w');
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
