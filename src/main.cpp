#include "geometry.h"
#include <iostream>

using namespace std;

int main() {
    Polygon subgrads = {
        Point(-4.0, -2.0),
        Point(1.0, -3.5),
        Point(4.5, 0.0),
        Point(2.5, 4.0),
        Point(-2.0, 3.5)
    };

    bool zero_in_hull = isZeroInConvexHull(subgrads);

    if (zero_in_hull) {
        cout << "ooptimum";
    } else {
        cout << "need moore steps";
    }

    // Polygon A = {
    //     Point(-3.5, 0.0),
    //     Point(-1.0, -2.5),
    //     Point(2.5, -1.5),
    //     Point(3.0, 2.0),
    //     Point(0.0, 3.5)
    // };

    // Polygon B = {
    //     Point(1.5, 4.0),
    //     Point(4.5, 0.5),
    //     Point(-2.0, -1.0)
    // };
    
    // Polygon Q = initQ0(A, B);
    
    // Point x = getCenter(Q);
  
    // SupportFunction rhoA(A), rhoB(B);

    // double best_F = numeric_limits<double>::infinity();

    // int iter = 0;
    // const int MAX_ITER = 1000;

    // while (iter < MAX_ITER) {
    //     FunctionValue fval = compute_F(x, A, B, rhoA, rhoB);
    //     vector<Point> subgrads = computeSubgradients(fval);

    //     if (fval.value < best_F) {
    //         best_F = fval.value;
    //     }

    //     cout << "Iter " << iter << ": x=(" << x.x << ", " << x.y 
    //          << "), F(x)=" << fval.value << ", best=" << best_F << endl;

    //     ExtremeDirections extreme = findExtremeDirections(subgrads);

    //     Q = clipPolygon(Q, x, extreme.min_angle);
    //     Q = clipPolygon(Q, x, extreme.max_angle);

    //     if (Q.empty()) {
    //         cout << "Error: Q became empty!" << endl;
    //         break;
    //     }

    //     for (size_t i = 0; i < Q.size(); ++i) {
    //         cout << "x:" << Q[i].x << " " << "y:" << Q[i].y << " ";

    //     }
    //     cout << endl;

    //     Point new_x = getCenter(Q);

    //     bool small_argument_change = (new_x - x).norm() < 1e-10;
    //     bool zero_in_hull = isZeroInConvexHull(subgrads);
        
    //     x = new_x;

    //     if (zero_in_hull || small_argument_change) {
    //         cout << "\nSTOP: ";
    //         if (zero_in_hull) cout << "0 in subgradient hull\n";
    //         if (small_argument_change) cout << "argument stopped changing\n";
    //         cout << "Best shift: (" << x.x << ", " << x.y << ")" << endl;
    //         break;
    //     }

    //     iter++;
    // }
    
    return 0;
}