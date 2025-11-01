#include <iostream>
#include <cmath>
#include <string>
#include <cassert>
#include <functional>
#include <set>

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

    double leftAreaAt0 = leftArea(alignedTrapezoid, 0.0);
    double leftAreaAt05 = leftArea(alignedTrapezoid, 0.5);
    double leftAreaAt1 = leftArea(alignedTrapezoid, 1.0);
    std::cout << "Left area at t=0: " << leftAreaAt0 << std::endl;
    std::cout << "Left area at t=0.5: " << leftAreaAt05 << std::endl;
    std::cout << "Left area at t=1: " << leftAreaAt1 << std::endl;

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
  // integral of max(0, min(x, x_r(y)) - x_l(y))
  //   when not truncated by vertical line at x
  static double fullSection(double y_a, double y_b, double h, double width1, double width2) {
    double dy = y_b - y_a;
    double area = width1*dy + (width2 - width1)*(y_b*y_b - y_a*y_a)/(2.0*h);
    return area;
  }

  // integral of max(0, min(x, x_r(y)) - x_l(y))
  //   when truncated by vertical line at x
  static double truncSection(double y_a, double y_b, double h, double x, double x0, double x1) {
    double dy = y_b - y_a;
    double area = (x - x0)*dy - (x1-x0)*(y_b*y_b - y_a*y_a)/(2.0*h);
    return area;
  }

  static double leftArea(const Trapezoid& trapezoid, double t) {
    assert(t >= 0 && t <= 1);
    auto a = trapezoid.getA();
    auto b = trapezoid.getB();
    auto c = trapezoid.getC();
    auto d = trapezoid.getD();
    double h = trapezoid.height();

    // find x at parameter t
    double x0 = std::min(a.x, d.x);
    double x1 = std::max(b.x, c.x);
    double x = x0 + t * (x1 - x0);

    // find vertical line intersections with left and right sides
    //   clamp within trapezoid
    double y_l = h*((x - a.x)/(d.x - a.x)); 
    y_l = std::clamp(y_l, 0.0, h);
    double y_r = h*((x - b.x)/(c.x - b.x));
    y_r = std::clamp(y_r, 0.0, h);

    // find unique points where area equation changes
    //   (completely on left, cut by vertical line, completely on right)
    std::set<double> breakpoints = {0.0, y_l, y_r, h};
    std::vector<double> bpVec(breakpoints.begin(), breakpoints.end());
    int intervals = breakpoints.size() - 1;

    // iterate over intervals where area equation is consistent
    //   and sum area contributions
    double area = 0.0;
    for (int i = 0; i < intervals; i++) {
      auto start = bpVec[i];
      auto end = bpVec[i+1];
      auto mid = 0.5 * (start + end);
      std::cout << "Interval [" << start << ", " << end << "], mid: " << mid << std::endl;

      double x_l_mid = a.x + (d.x - a.x) / h * mid;
      double x_r_mid = b.x + (c.x - b.x) / h * mid;
      std::cout << "  x_l_mid: " << x_l_mid << ", x_r_mid: " << x_r_mid << std::endl;
      if (x <= x_l_mid) {
        std::cout << "0 area" << std::endl;
        continue;
      } else if (x >= x_r_mid) {
        std::cout << "full area" << std::endl;
        area += fullSection(start, end, h,
            Point::distance(a, b), Point::distance(d, c));
      } else {
        std::cout << "trunc area" << std::endl;
        area += truncSection(start, end, h, x,
            a.x, d.x);
      }
    }

    return area;
  }

  static double rootFunc(const Trapezoid& trapezoid, double t , double targetArea) {
    return leftArea(trapezoid, t) - targetArea;
  }

  static std::pair<Point, Point> calculateEF(const Trapezoid& trapezoid, double t) {
    assert(t >= 0 && t <= 1);
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

    double x0_l = std::min(a.x, d.x);
    double x1_l = std::max(a.x, d.x);
    double x0_r = std::min(b.x, c.x);
    double x1_r = std::max(b.x, c.x);

    // vertical line intersects left side
    if (x > x0_l && x < x1_l) {
      double y_l = h*((x - a.x)/(d.x - a.x)); 
      Point l = a.x < d.x ? a : d;
      Point r = a.x < d.x ? d : a;
      // determine which is top and which is bottom
      //   depending on the slope sign
      if (l.y - r.y < 0) {
        y_top = y_l;
      } else {
        y_bottom = y_l;
      }
    }

    // vertical line intersects right side
    if (x > x0_r && x < x1_r) {
      double y_r = h*((x - b.x)/(c.x - b.x));
      Point l = b.x < c.x ? b : c;
      Point r = b.x < c.x ? c : b;
      if (l.y - r.y < 0) {
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
