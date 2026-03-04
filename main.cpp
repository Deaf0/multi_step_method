#include "geometry.h"
#include <iostream>

using namespace std;

int main() {
    Polygon A = {
        Point(-2, -1),
        Point(2, -1),
        Point(2, 1),
        Point(-2, 1)
    };

// Прямоугольник B (подвижный) - повернут на 90 градусов по сути
    Polygon B = {
        Point(-0.5, -1.5),
        Point(0.5, -1.5),
        Point(0.5, 1.5),
        Point(-0.5, 1.5)
    }; 
    
    Polygon Q = initQ0(A, B);
    
    Point x = getCenter(Q);
  
    SupportFunction rhoA(A), rhoB(B);

    double best_F = numeric_limits<double>::infinity();

    int iter = 0;
    const int MAX_ITER = 1000;

    while (iter < MAX_ITER) {
        FunctionValue fval = compute_F(x, A, B, rhoA, rhoB);
        vector<Point> subgrads = computeSubgradients(fval);

        if (fval.value < best_F) {
            best_F = fval.value;
        }

        cout << "Iter " << iter << ": x=(" << x.x << ", " << x.y 
             << "), F(x)=" << fval.value << ", best=" << best_F << endl;

        ExtremeDirections extreme = findExtremeDirections(subgrads);

        Q = clipPolygon(Q, x, extreme.min_angle);
        Q = clipPolygon(Q, x, extreme.max_angle);

        if (Q.empty()) {
            cout << "Error: Q became empty!" << endl;
            break;
        }

        Point new_x = getCenter(Q);

        bool small_argument_change = (new_x - x).norm() < 1e-10;
        bool zero_in_hull = isZeroInConvexHull(subgrads);
        
        x = new_x;

        if (zero_in_hull || small_argument_change) {
            cout << "\nSTOP: ";
            if (zero_in_hull) cout << "0 in subgradient hull\n";
            if (small_argument_change) cout << "argument stopped changing\n";
            cout << "Best shift: (" << x.x << ", " << x.y << ")" << endl;
            break;
        }

        iter++;
    }
    
    return 0;
}