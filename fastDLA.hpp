#include <vector>

class ClusterTree
{
public:
  ClusterTree(const double maxRadius, const double minLength = 4) : root_(nullptr)
  {
    minSideLength = minLength;
    maxDepth_ = 1 + ceil_log2(maxRadius / minSideLength);

    Point point;
    point.x=0.0;
    point.y=0.0;
    insertPoint(point);

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

  struct Point
  {
    double x;
    double y;
  };
  
  struct Node
  {
    Node* NW;
    Node* NE;
    Node* SW;
    Node* SE;
    std::vector<Point> points;
  };

  Node* root_;
  const int maxDepth;
  const double minSideLength;
  double startDist;
}
