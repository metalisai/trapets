#include <iostream>
#include <cmath>
#include <string>
#include <cassert>
#include <functional>

struct Point {
  double x;
  double y;

  static double distance(const Point& a, const Point& b) {
    return std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y));
  }

  std::string toString() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
  }
};

class Trapezoid {
public:
  Trapezoid(const Point& a, const Point& b, const Point& c, const Point& d)
      : a(a), b(b), c(c), d(d) {}

  // validate if the points form a valid trapezoid
  static bool validate(const Point& a, const Point& b, const Point& c, const Point& d) {
    // check for degenerate vertical lines
    if (b.x == a.x || d.x == c.x) {
      return false;
    }
    // check if the top and bottom sides are parallel
    double slope1 = (b.y - a.y) / (b.x - a.x);
    double slope2 = (d.y - c.y) / (d.x - c.x);
    return std::abs(slope1 - slope2) < 1e-9;
  }

  double area() const {
    double base1 = Point::distance(a, b);
    double base2 = Point::distance(c, d);
    double height = std::abs(a.y - c.y);
    return 0.5 * (base1 + base2) * height;
  }

  std::string toString() const {
    return "Trapezoid(" + a.toString() + ", " + b.toString() + ", " +
           c.toString() + ", " + d.toString() + ")";
  }
private:
  Point a, b, c, d;
};

class BisectionRootFinder {
public:
  static double solve(std::function<double(double)> f,
      double a, double b, double tol = 1e-9, int maxIter = 100) {
    throw std::runtime_error("Not implemented");
  }
};

class Challenge {
  public:
  static std::pair<Point, Point> calculateEF(const Trapezoid& trapezoid) {
    auto area = trapezoid.area();
    auto halfArea = area / 2.0;

    // TODO: check if bases are horizontal
    //   if not, rotate the trapezoid to make them horizontal
    Trapezoid alignedTrapezoid = trapezoid;

    double t = BisectionRootFinder::solve(
      [&](double t) {
        return rootFunc(alignedTrapezoid, t, halfArea);
      },
      0.0, 1.0
    );

    return {{0, 0}, {0, 0}}; // placeholder implementation
  }
  protected:
  static double leftArea(const Trapezoid& trapezoid, double t) {
    assert(t >= 0 && t <= 1);
    return 0.0; // placeholder implementation
  }

  static double rootFunc(const Trapezoid& trapezoid, double t , double targetArea) {
    return leftArea(trapezoid, t) - targetArea;
  }
};

int main(void) {
  Point a, b, c, d;
  std::cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;
  bool valid = Trapezoid::validate(a, b, c, d);
  if (!valid) {
    std::cerr << "Invalid trapezoid points ABCD. AB and CD must be parallel and non-zero length." << std::endl;
    return 1;
  }
  Trapezoid trapezoid(a, b, c, d);
  std::cout << trapezoid.toString() << std::endl;
  std::cout << "Area: " << trapezoid.area() << std::endl;

  std::cout << "Exiting" << std::endl;
  return 0;
}
