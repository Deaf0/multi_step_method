#include "geometry.h"
#include <stdexcept>


// === Point ===

Point::Point(double x, double y) : x(x), y(y) {}

Point Point::operator+(const Point& other) const {
    return Point(x + other.x, y + other.y);
}

Point Point::operator-(const Point& other) const {
    return Point(x - other.x, y - other.y);
}

Point Point::operator*(double scalar) const {
    return Point(x * scalar, y * scalar);
}

double Point::dot(const Point& other) const {
    return x * other.x + y * other.y;
}

double Point::norm() const {
    return sqrt(x*x + y*y);
}

// === SupportFunction ===

SupportFunction::SupportFunction(const Polygon& p) : poly(p) {}

double SupportFunction::operator()(const Point& l) const {
    double max_val = -numeric_limits<double>::infinity();
    for (const auto& v : poly) {
        max_val = max(max_val, l.dot(v));
    }
    return max_val;
}

// === BoundingBox ===

BoundingBox::BoundingBox(const Polygon& poly) {
    min.x = min.y = numeric_limits<double>::infinity();
    max.x = max.y = -numeric_limits<double>::infinity();

    for (const auto& p : poly) {
        min.x = std::min(min.x, p.x);
        min.y = std::min(min.y, p.y);
        max.x = std::max(max.x, p.x);
        max.y = std::max(max.y, p.y);
    }
}

// === FunctionValue ===

FunctionValue::FunctionValue(double v) : value(v) {}

// === ActiveVertices ===

ActiveVertices findActiveVertices(const Point& l, const Polygon& A, const Polygon& B) {
    ActiveVertices result;
    result.rhoA_val = -numeric_limits<double>::infinity();
    result.rhoB_val = -numeric_limits<double>::infinity();

    for (size_t i = 0; i < A.size(); ++i) {
        double val = l.dot(A[i]);
        if (val > result.rhoA_val) {
            result.rhoA_val = val;
            result.a_idx = i;
        }
    }

    for (size_t j = 0; j < B.size(); ++j) {
        double val = l.dot(B[j]);
        if (val > result.rhoB_val) {
            result.rhoB_val = val;
            result.b_idx = j;
        }
    }

    return result;
}

// === Phi functions ===

double computePhiI(const Point& x, const Point& a_i, const Point& l, const SupportFunction& rhoB) {
    return l.dot(a_i - x) - rhoB(l);
}

double computePhiJ(const Point& x, const Point& b_j, const Point& l, const SupportFunction& rhoA) {
    return l.dot(b_j + x) - rhoA(l);  
}

// === compute_F ===

FunctionValue compute_F(const Point& x,
                        const Polygon& A,
                        const Polygon& B,
                        const SupportFunction& rhoA,
                        const SupportFunction& rhoB) {

    FunctionValue result(-numeric_limits<double>::infinity());

    for (size_t i = 0; i < A.size(); ++i) {
        Point vec = A[i] - x;
        double len = vec.norm();
        if (len < 1e-10) continue;

        Point l = vec * (1.0 / len);
        double phi = computePhiI(x, A[i], l, rhoB);

        if (phi > result.value + 1e-10) {
            result.value = phi;
            result.directions.clear();
            result.directions.push_back(l);
            result.active_i.clear();
            result.active_j.clear();
            result.active_i.push_back(i);
            auto active = findActiveVertices(l, A, B);
            result.active_j.push_back(active.b_idx);
        } else if (std::abs(phi - result.value) < 1e-10) {
            result.directions.push_back(l);
            result.active_i.push_back(i);
            auto active = findActiveVertices(l, A, B);
            result.active_j.push_back(active.b_idx);
        }
    }

    for (size_t j = 0; j < B.size(); ++j) {
        Point vec = B[j] + x;
        double len = vec.norm();
        if (len < 1e-10) continue;

        Point l = vec * (1.0 / len);
        double phi = computePhiJ(x, B[j], l, rhoA);

        if (phi > result.value + 1e-10) {
            result.value = phi;
            result.directions.clear();
            result.directions.push_back(l);
            result.active_i.clear();
            result.active_j.clear();
            result.active_j.push_back(j);
            auto active = findActiveVertices(l, A, B);
            result.active_i.push_back(active.a_idx);
        } else if (std::abs(phi - result.value) < 1e-10) {
            result.directions.push_back(l);
            result.active_j.push_back(j);
            auto active = findActiveVertices(l, A, B);
            result.active_i.push_back(active.a_idx);
        }
    }

    return result;
}

// === computeSubgradients ===

vector<Point> computeSubgradients(const FunctionValue& fval) {
    vector<Point> subgrads;

    for (size_t idx = 0; idx < fval.directions.size(); ++idx) {
        const Point& l = fval.directions[idx];

        if (idx < fval.active_i.size() && fval.active_i[idx] >= 0) {
            subgrads.push_back(l * (-1.0));
        }

        if (idx < fval.active_j.size() && fval.active_j[idx] >= 0) {
            subgrads.push_back(l);
        }
    }

    return subgrads;
}

// === initQ0 ===

Polygon initQ0(const Polygon& A, const Polygon& B) {
    BoundingBox box_a(A), box_b(B);

    Point bottom_left(
        box_b.max.x - box_a.min.x,
        box_b.max.y - box_a.min.y
    );

    Point top_right(
        box_a.max.x - box_b.min.x,
        box_a.max.y - box_b.min.y
    );

    Polygon Q0;
    Q0.push_back(bottom_left);
    Q0.push_back(Point(top_right.x, bottom_left.y));
    Q0.push_back(top_right);
    Q0.push_back(Point(bottom_left.x, top_right.y));

    return Q0;
}

// === getCenter ===

Point getCenter(const Polygon& poly) {
    Point center(0, 0);
    for (const auto& p : poly) {
        center = center + p;
    }
    return center * (1.0 / poly.size());
}

// === signedDistance ===

double signedDistance(const Point& p, const Point& n, const Point& point) {
    return n.dot(point - p);
}

// === clipPolygon ===

Polygon clipPolygon(const Polygon& poly, const Point& p, const Point& n) {
    Polygon result;

    if (poly.empty()) return result;

    size_t size = poly.size();
    for (size_t i = 0; i < size; ++i) {
        const Point& curr = poly[i];
        const Point& next = poly[(i + 1) % size];

        double curr_dist = signedDistance(p, n, curr);
        double next_dist = signedDistance(p, n, next);

        if (curr_dist <= 0) {
            result.push_back(curr);
        }

        if (curr_dist * next_dist < 0) {
            double t = curr_dist / (curr_dist - next_dist);
            Point intersection(
                curr.x + t * (next.x - curr.x),
                curr.y + t * (next.y - curr.y)
            );
            result.push_back(intersection);  
        }
    }

    return result;
}

// === findExtremeDirections ===

ExtremeDirections findExtremeDirections(const vector<Point>& vectors) {
    if (vectors.empty()) {
        throw runtime_error("no subgradients");
    }

    double min_angle = atan2(vectors[0].y, vectors[0].x);
    double max_angle = min_angle;
    Point min_vec = vectors[0];
    Point max_vec = vectors[0];

    for (const auto& v : vectors) {
        double angle = atan2(v.y, v.x);
        if (angle < min_angle) {
            min_angle = angle;
            min_vec = v;
        }
        if (angle > max_angle) {
            max_angle = angle;
            max_vec = v;
        }
    }

    return {min_vec, max_vec};
}

// === isZeroInConvexHull (временна€ заглушка) ===

bool isZeroInConvexHull(const vector<Point>& vectors) {
    // TODO: реализовать GJK или метод углов
    // ѕока всегда возвращаем false, чтобы алгоритм не останавливалс€ раньше времени
    return false;
}