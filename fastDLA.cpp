#include "treeDLA.hpp"
#include "gridDLA.hpp"
#include <iostream>

int main(int argc, char** argv)
{
  int n, seed;
  char filename[1000];
  char method=0;

  seed=time(NULL);
  filename[0]='\0';

  switch(argc-1)
    {
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
	ClusterGrid grid(n);

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
		grid.writePoints(fp);
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
