#include "treeDLA.hpp"
#include "gridDLA.hpp"
#include <iostream>
#include <stdlib.h>

int compareDouble (const void * a, const void * b)
{
  if ( *(double*)a <  *(double*)b ) return -1;
  if ( *(double*)a >  *(double*)b ) return 1;
  if ( *(double*)a == *(double*)b ) return 0;
}

int main(int argc, char** argv)
{
  int n, seed;
  char filename[1000];
  char method=0;
  char resultType='w';
  double noiseReductionFactor=1;

  int distsCount=10000;
  double dists[10000];//must match above line

  seed=time(NULL);
  filename[0]='\0';

  switch(argc-1)
    {
    case 6:
      sscanf(argv[6],"%lf",&noiseReductionFactor);
    case 5:
      sscanf(argv[5],"%c",&resultType);
    case 4:
      sscanf(argv[4],"%d",&seed);
    case 3:
      sscanf(argv[3],"%s",&filename);
    case 2:
      sscanf(argv[2],"%c",&method);
    case 1:
      sscanf(argv[1],"%d",&n);
      break;
    default:
      fprintf(stderr,"First 2 arguments (Cluster Size and method) are mandatory.\nThird argument (Filename to save to) and fourth argument (seed) are optional.\n");
      return 1;
    }

  switch (method)
    {
    case 'g':
      {
	ClusterGrid grid(n,seed,noiseReductionFactor,6+noiseReductionFactor*18,1+noiseReductionFactor);

	if (resultType=='d')
	  {
	    for (int i=0;i<distsCount;i++)
	      {
		dists[i]=grid.getNormalisedDist();
	      }
	    qsort(dists,distsCount,sizeof(double),compareDouble);
	  }

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
		if (resultType=='w')
		  {
		    grid.writePoints(fp);
		  }
		else if (resultType=='d')
		  {
		    fwrite(dists,sizeof(double),distsCount,fp);
		  }
		fclose(fp);
		printf("Written cluster to file\n");
	      }
	  }
	else
	  {
	    printf("Not saving points as no filename given\n");
	  }
	break;
      }
    case 't':
      {
	ClusterTree tree(n);

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
		tree.writePoints(fp);
		fclose(fp);
		printf("Written cluster to file\n");
	      }
	  }
	else
	  {
	    printf("Not saving points as no filename given\n");
	  }
	break;
      }
    default:
      {
	printf("method (second argument) must be either g for grid or t for tree\n");
	return 1;
      }
    }



  return 0;
}
