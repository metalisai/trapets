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

  double height() const {
    // TODO: handle non-horizontal bases
    return std::abs(a.y - c.y);
  }

  double area() const {
    double base1 = Point::distance(a, b);
    double base2 = Point::distance(c, d);
    double height = this->height();
    return 0.5 * (base1 + base2) * height;
  }

  std::string toString() const {
    return "Trapezoid(" + a.toString() + ", " + b.toString() + ", " +
           c.toString() + ", " + d.toString() + ")";
  }

  const Point& getA() const { return a; }
  const Point& getB() const { return b; }
  const Point& getC() const { return c; }
  const Point& getD() const { return d; }

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

    /*double t = BisectionRootFinder::solve(
      [&](double t) {
        return rootFunc(alignedTrapezoid, t, halfArea);
      },
      0.0, 1.0
    );*/
    double t = 0.5;
    auto ef = calculateEF(alignedTrapezoid, t);

    return ef;
  }
  protected:
  static double leftArea(const Trapezoid& trapezoid, double t) {
    assert(t >= 0 && t <= 1);
    return 0.0; // placeholder implementation
  }

  static double rootFunc(const Trapezoid& trapezoid, double t , double targetArea) {
    return leftArea(trapezoid, t) - targetArea;
  }

  static std::pair<Point, Point> calculateEF(const Trapezoid& trapezoid, double t) {
    auto a = trapezoid.getA();
    auto b = trapezoid.getB();
    auto c = trapezoid.getC();
    auto d = trapezoid.getD();

    double h = trapezoid.height();
    double x0 = std::min(a.x, d.x);
    double x1 = std::max(b.x, c.x);
    double x = x0 + t * (x1 - x0);
    std::cout << "x: " << x << std::endl;

    double y_bottom = a.y;
    double y_top = y_bottom + h;

    double slope_l = ((x - a.x)/(d.x - a.x));
    double slope_r = ((x - b.x)/(c.x - b.x));

    double x0_l = std::min(a.x, d.x);
    double x1_l = std::max(a.x, d.x);
    double x0_r = std::min(b.x, c.x);
    double x1_r = std::max(b.x, c.x);

    double y_l = h*slope_l; 
    double y_r = h*slope_r;

    std::cout << "y_l: " << y_l << ", y_r: " << y_r << std::endl;

    if (x > x0_l && x < x1_l) {
      if (slope_l < 0) {
        y_bottom = y_l;
      } else {
        y_top = y_l;
      }
    }

    if (x > x0_r && x < x1_r) {
      if (slope_r >= 0) {
        y_bottom = y_r;
      } else {
        y_top = y_r;
      }
    }

    Point e = {x, y_bottom};
    Point f = {x, y_top};
    return {e, f};
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
  auto [e, f] = Challenge::calculateEF(trapezoid);
  std::cout << "E: " << e.toString() << std::endl;
  std::cout << "F: " << f.toString() << std::endl;

  std::cout << "Exiting" << std::endl;
  return 0;
}
