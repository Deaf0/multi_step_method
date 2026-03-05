#pragma once

#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

using namespace std;

struct Point {
    double x, y;

    Point(double x = 0, double y = 0);

    Point operator+(const Point& other) const;
    Point operator-(const Point& other) const;
    Point operator*(double scalar) const;

    double dot(const Point& other) const;
    double norm() const;
};

using Polygon = vector<Point>;

struct SupportFunction {
    const Polygon& poly;

    SupportFunction(const Polygon& p);
    double operator()(const Point& l) const;
};

struct BoundingBox {
    Point min, max;

    BoundingBox(const Polygon& poly);
};

struct FunctionValue {
    double value;
    vector<Point> directions;
    vector<int> active_i;
    vector<int> active_j;

    FunctionValue(double v = 0);
};

struct ActiveVertices {
    int a_idx;
    int b_idx;
    double rhoA_val;
    double rhoB_val;
};

struct ExtremeDirections {
    Point min_angle;
    Point max_angle;
};


ActiveVertices findActiveVertices(const Point& l, const Polygon& A, const Polygon& B);

double computePhiI(const Point& x, const Point& a_i, const Point& l, const SupportFunction& rhoB);
double computePhiJ(const Point& x, const Point& b_j, const Point& l, const SupportFunction& rhoA);

FunctionValue compute_F(const Point& x,
                        const Polygon& A,
                        const Polygon& B,
                        const SupportFunction& rhoA,
                        const SupportFunction& rhoB);

vector<Point> computeSubgradients(const FunctionValue& fval);

Polygon initQ0(const Polygon& A, const Polygon& B);
Point getCenter(const Polygon& poly);

double signedDistance(const Point& p, const Point& n, const Point& point);
Polygon clipPolygon(const Polygon& poly, const Point& p, const Point& n);

ExtremeDirections findExtremeDirections(const vector<Point>& vectors);
bool isZeroInConvexHull(const vector<Point>& vectors);  
