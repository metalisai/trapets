#include <iostream>
#include <cmath>
#include <string>
#include <cassert>
#include <functional>
#include <set>
#include <vector>
#include <algorithm>

#include <fstream>

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

struct Mat2x2 {
  double m00, m01;
  double m10, m11;
};

Point operator*(const Mat2x2& mat, const Point& p) {
  return {
    mat.m00 * p.x + mat.m01 * p.y,
    mat.m10 * p.x + mat.m11 * p.y
  };
}

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
    return std::abs(slope1 - slope2) < 1e-4;
  }

  // return a rotated copy of the trapezoid
  Trapezoid rotated(double angleRad) const {
    double cosA = std::cos(angleRad);
    double sinA = std::sin(angleRad);
    Mat2x2 rotationMatrix = {cosA, -sinA, sinA, cosA};
    return Trapezoid(
      rotationMatrix * a,
      rotationMatrix * b,
      rotationMatrix * c,
      rotationMatrix * d
    );
  }

  // compute height of trapezoid as distance between parallel sides AB and CD
  double height() const {
    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double length = std::sqrt(dx*dx + dy*dy);
    double height = std::abs((d.x - a.x)*dy - (d.y - a.y)*dx) / length;
    return height;
  }

  // compute area of trapezoid
  double area() const {
    double base1 = Point::distance(a, b);
    double base2 = Point::distance(c, d);
    double height = this->height();
    return 0.5 * (base1 + base2) * height;
  }

  // return string representation
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

// simple bisection root finder
class BisectionRootFinder {
  public:

  static double solve(std::function<double(double)> f,
      double a, double b, double tol = 1e-9, int maxIter = 100) {
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0) {
      throw std::invalid_argument("Function must have opposite signs at the endpoints.");
    }
    for (int iter = 0; iter < maxIter; iter++) {
      double c = 0.5 * (a + b);
      double fc = f(c);
      if (std::abs(fc) < tol || (b - a) / 2 < tol) {
        return c;
      }
      if (fa * fc < 0) {
        b = c;
        fb = fc;
      } else {
        a = c;
        fa = fc;
      }
    }
    throw std::runtime_error("Maximum iterations reached without convergence.");
  }
};

class Challenge {
  public:

  // solve the challenge: find points E and F that split trapezoid area in half
  static std::pair<Point, Point> calculateEF(const Trapezoid& trapezoid) {

    double angleRad = std::atan2(trapezoid.getB().y - trapezoid.getA().y,
                                 trapezoid.getB().x - trapezoid.getA().x);
    angleRad = std::abs(angleRad) < 1e-9 ? 0.0 : angleRad;
    Trapezoid alignedTrapezoid = trapezoid.rotated(-angleRad);

    auto area = alignedTrapezoid.area();
    auto halfArea = area / 2.0;

    double t = BisectionRootFinder::solve(
      [&](double t) {
        return rootFunc(alignedTrapezoid, t, halfArea);
      },
      0.0, 1.0
    );

    auto ef = calculateEFAtT(alignedTrapezoid, t);

    // rotate points E and F back to original orientation
    double cosA = std::cos(angleRad);
    double sinA = std::sin(angleRad);
    Mat2x2 rot = {
      cosA, -sinA,
      sinA,  cosA
    };
    ef.first = rot * ef.first;
    ef.second = rot * ef.second;

    return ef;
  }
  private:

  // integral of max(0, min(x, x_r(y)) - x_l(y)) over vertical section
  //   when not truncated by vertical line at x
  static double fullSection(double y_a, double y_b, double h, double width1, double width2) {
    double dy = y_b - y_a;
    double area = width1*dy + (width2 - width1)*(y_b*y_b - y_a*y_a)/(2.0*h);
    return area;
  }

  // integral of max(0, min(x, x_r(y)) - x_l(y)) over vertical section
  //   when truncated by vertical line at x
  static double truncSection(double y_a, double y_b, double h, double x, double x0, double x1) {
    double dy = y_b - y_a;
    double area = (x - x0)*dy - (x1-x0)*(y_b*y_b - y_a*y_a)/(2.0*h);
    return area;
  }

  // compute area to the left of vertical line at parameter t
  //   trapezoid is assumed to have horizontal bases
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
    assert(intervals >= 1);

    // iterate over intervals where area equation is consistent
    //   and sum area
    double area = 0.0;
    for (int i = 0; i < intervals; i++) {
      auto start = bpVec[i];
      auto end = bpVec[i+1];
      auto mid = 0.5 * (start + end);

      // left and right x positions at mid height of interval
      double x_l_mid = a.x + (d.x - a.x) / h * mid;
      double x_r_mid = b.x + (c.x - b.x) / h * mid;
      if (x <= x_l_mid) {
        continue;
      } else if (x >= x_r_mid) {
        area += fullSection(start, end, h,
            Point::distance(a, b), Point::distance(d, c));
      } else {
        area += truncSection(start, end, h, x,
            a.x, d.x);
      }
    }

    return area;
  }

  // function for root finding
  static double rootFunc(const Trapezoid& trapezoid, double t, double targetArea) {
    return leftArea(trapezoid, t) - targetArea;
  }

  // calculate points E and F at specified horizontal position
  //   t is in [0, 1]
  static std::pair<Point, Point> calculateEFAtT(const Trapezoid& trapezoid, double t) {
    assert(t >= 0 && t <= 1);
    auto a = trapezoid.getA();
    auto b = trapezoid.getB();
    auto c = trapezoid.getC();
    auto d = trapezoid.getD();

    double h = trapezoid.height();
    double x0 = std::min(a.x, d.x);
    double x1 = std::max(b.x, c.x);
    double x = x0 + t * (x1 - x0);

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
  // read input
  Point a, b, c, d;
  std::cout << "Enter trapezoid points A B C D as x1 y1 x2 y2 x3 y3 x4 y4:" << std::endl;
  std::cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;

  // validate input
  bool valid = Trapezoid::validate(a, b, c, d);
  if (!valid) {
    std::cerr << "Invalid trapezoid points ABCD. AB and CD must be parallel and non-zero length." << std::endl;
    return 1;
  }

  // solve challenge
  Trapezoid trapezoid(a, b, c, d);
  std::cout << trapezoid.toString() << std::endl;
  std::cout << "Full trapezoid area: " << trapezoid.area() << std::endl;
  try {
    auto [e, f] = Challenge::calculateEF(trapezoid);

    // output results to stdout
    std::cout << "E: " << e.toString() << std::endl;
    std::cout << "F: " << f.toString() << std::endl;

    // write output to file
    //   first line - input trapezoid points
    //   second line - output points E and F
    std::cout << "Writing file trapezoid.txt" << std::endl;
    std::ofstream ofs("trapezoid.txt");
    ofs << a.x << " " << a.y << " "
        << b.x << " " << b.y << " "
        << c.x << " " << c.y << " "
        << d.x << " " << d.y << std::endl
        << e.x << " " << e.y << " "
        << f.x << " " << f.y << std::endl;
    ofs.close();
  } catch (const std::exception& ex) {
    std::cerr << "Error solving challenge: " << ex.what() << std::endl;
    return 1;
  }

  std::cout << "Exiting" << std::endl;
  return 0;
}
