//
// Created by nack2 on 04.06.2025.
//
#include "../include/nonLinearEquation.h"

#include <algorithm>
#include <iostream>

double metoda_newton(double (*f)(double), double y0){
    double y1 = y0;
    for (int i = 0; i < ITERATION_LIMIT; ++i){
        count_newton++;
        double f_y0 = f(y0);
        double prim_y0 = prim(f, y0);

        // std::cout<< "y0: "<<y0<<"\n";
        // std::cout<<y0<<"\n";
        // std::cout<< "f(y0): "<<f(y0)<<"\n";

        if (std::isnan(f_y0) || std::isinf(f_y0) || std::isnan(prim_y0) || std::isinf(prim_y0)) {
            // std::cerr << "Nieprawidlowa wartosc: f(y0)=" << f_y0 << ", prim=" << prim_y0 << "\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (std::abs(prim_y0) < EPSILON) {
            // std::cerr << "WARNING: Pochodna bliska zeru: prim(y0)=" << prim_y0 << "\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        y1 = y0 - f_y0 / prim_y0;

        if (std::abs(f(y1)) < EPSILON || std::abs(y1 - y0) < EPSILON) {
            return y1;
        }

        y0 = y1;
    }

    std::cerr << "WARNING: Out of iteration limit\n";
    // return std::numeric_limits<double>::quiet_NaN();
    return y0;
}

double metoda_newton(double (*f)(double), double (*f_prim)(double), double y0){
    double y1 = y0;
    for (int i = 0; i < ITERATION_LIMIT; ++i){
        double f_y0 = f(y0);
        double prim_y0 = f_prim(y0);

        if (std::isnan(f_y0) || std::isnan(prim_y0) || std::isinf(prim_y0)) {
            // std::cerr << "Nieprawidlowa wartosc: f(y0)=" << f_y0 << ", prim=" << prim_y0 << "\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (std::abs(prim_y0) < DERIVATIVE_EPSILON) {
            // prim_y0 = std::copysign(DERIVATIVE_EPSILON, prim_y0);
            std::cerr << "Pochodna zbyt bliska 0.\n";
            return std::numeric_limits<double>::quiet_NaN();
        }

        double step = -f_y0 / prim_y0;
        if (std::abs(step) > MAX_STEP)
            step = std::copysign(MAX_STEP, step);

        y1 = y0 + step;

        if (std::abs(y1 - y0) <= EPSILON){
            return y1;
        }

        y0 = y1;
    }

    std::cerr << "Przekroczono limit\n";
    // return std::numeric_limits<double>::quiet_NaN();
    return y0;
}
double metoda_bisekcji(std::function<double(double)> f, double a, double b, std::vector<double>& iteracje) {
    double fa_val = f(a); // Zmieniono nazwę, żeby nie kolidowała z parametrem (choć tu nie ma problemu)
    double fb_val = f(b);

    if (std::fabs(fa_val) < EPSILON) {
        iteracje.push_back(a);
        return a;
    }
    if (std::fabs(fb_val) < EPSILON) {
        iteracje.push_back(b);
        return b;
    }

    if (!std::isfinite(fa_val) || !std::isfinite(fb_val) || fa_val * fb_val > 0)
        return NAN;

    double current_a = a; // Używamy lokalnych kopii do modyfikacji
    double current_b = b;
    double fa_current = fa_val; // fa_current odpowiada current_a

    double mid = current_a;
    for (int i = 0; i < ITERATION_LIMIT; ++i) {
        mid = (current_a + current_b) / 2.0;
        // Zabezpieczenie przed sytuacją, gdy mid przestaje się zmieniać z powodu precyzji
        if (mid == current_a || mid == current_b) {
            break;
        }

        double fmid = f(mid);
        iteracje.push_back(mid);

        if (std::fabs(fmid) < EPSILON || std::fabs(current_b - current_a) / 2.0 < EPSILON) // Użyto (b-a)/2 dla błędu przybliżenia
            break;

        if (fa_current * fmid < 0) {
            current_b = mid;
            // fb_val = fmid; // Nie jest potrzebne, bo fb_val nie jest używane w kolejnej iteracji w tej gałęzi
        } else {
            current_a = mid;
            fa_current = fmid;
        }
    }
    return mid;
}

double regulaFalsi(function<double(double)> f, double a, double b, vector<double>& iteracje) {
    double fa_val = f(a);
    double fb_val = f(b);

    if (std::fabs(fa_val) < EPSILON) {
        iteracje.push_back(a);
        return a;
    }
    if (std::fabs(fb_val) < EPSILON) {
        iteracje.push_back(b);
        return b;
    }

    if (!isfinite(fa_val) || !isfinite(fb_val) || fa_val * fb_val > 0)
        return NAN;

    double current_a = a;
    double current_b = b;
    double fa_current = fa_val;
    double fb_current = fb_val;

    double x = current_a;
    double prev_x = current_a + 2 * EPSILON; // Inicjalizacja różna od x

    for (int i = 0; i < ITERATION_LIMIT; ++i) {
        if (std::fabs(fb_current - fa_current) < std::numeric_limits<double>::epsilon()) {
            break;
        }
        prev_x = x;
        x = (current_a * fb_current - current_b * fa_current) / (fb_current - fa_current);
        if (x == prev_x) { // x nie zmienia się
            if (std::fabs(f(x)) < EPSILON) break;
            break;
        }
        double fx = f(x);
        iteracje.push_back(x);
        if (std::fabs(fx) < EPSILON || std::fabs(x - prev_x) < EPSILON)
            break;
        if (fa_current * fx < 0) {
            current_b = x;
            fb_current = fx;
        } else {
            current_a = x;
            fa_current = fx;
        }
    }
    return x;
}

double metoda_siecznych(std::function<double(double)> f, double x0, double x1, std::vector<double>& iteracje) {
    double fx0 = f(x0);
    double fx1 = f(x1);

    if (std::fabs(fx0) < EPSILON) {
        iteracje.push_back(x0);
        return x0;
    }
    if (std::fabs(fx1) < EPSILON) {
        iteracje.push_back(x1);
        return x1;
    }

    if (!std::isfinite(fx0) || !std::isfinite(fx1)) {
        return NAN;
    }

    double current_x0 = x0;
    double current_x1 = x1;
    double current_fx0 = fx0;
    double current_fx1 = fx1;

    double x_next = current_x1;

    for (int i = 0; i < ITERATION_LIMIT; ++i) {
        if (std::fabs(current_fx1 - current_fx0) < std::numeric_limits<double>::epsilon()) {
            break;
        }

        x_next = current_x1 - current_fx1 * (current_x1 - current_x0) / (current_fx1 - current_fx0);

        // Jeśli x_next jest takie samo jak current_x1 (np. z powodu precyzji)
        if (x_next == current_x1) {
             // Sprawdź, czy to dlatego, że jesteśmy blisko rozwiązania
            if (std::fabs(f(x_next)) < EPSILON) break;
            // Inaczej metoda utknęła
            break;
        }

        iteracje.push_back(x_next);
        double fx_next = f(x_next);

        if (!std::isfinite(fx_next)) {
            return NAN;
        }

        if (std::fabs(fx_next) < EPSILON || std::fabs(x_next - current_x1) < EPSILON) {
            return x_next;
        }

        current_x0 = current_x1;
        current_fx0 = current_fx1;
        current_x1 = x_next;
        current_fx1 = fx_next;
    }
    return x_next;
}


void findAllRoots(function<double(double)> f, double start, double end, const string& fname) {
    cout << "\n--- " << fname << " ---\n";
    cout << fixed << setprecision(10);

    vector<double> found_roots_falsi;
    vector<double> found_roots_bisekcja;
    vector<double> found_roots_sieczne;

    // Metoda Regula Falsi
    for (double a_loop = start; a_loop < end; a_loop += STEP) {
        double b_loop = a_loop + STEP;
        if (b_loop > end) b_loop = end;

        double fa = f(a_loop), fb = f(b_loop);
        if (!isfinite(fa) || !isfinite(fb)) continue;

        if (fa * fb <= 0) {
            vector<double> iteracje;
            double root = regulaFalsi(f, a_loop, b_loop, iteracje);
            if (!isnan(root) && root >= a_loop - EPSILON && root <= b_loop + EPSILON) { // Sprawdź czy pierwiastek jest w sensownym zakresie
                found_roots_falsi.push_back(root);
            }
        }
         // Sprawdzenie, czy sam 'a_loop' nie jest pierwiastkiem (jeśli fa*fb > 0)
        else if (std::fabs(fa) < EPSILON) {
             found_roots_falsi.push_back(a_loop);
        }
        if (a_loop + STEP >= end && std::fabs(fb) < EPSILON) { // Dla ostatniego punktu 'b_loop'
             found_roots_falsi.push_back(b_loop);
        }
    }
    ranges::sort(found_roots_falsi.begin(), found_roots_falsi.end());
    found_roots_falsi.erase(unique(found_roots_falsi.begin(), found_roots_falsi.end(),
                                  [](double r1, double r2){ return std::fabs(r1-r2) < EPSILON*10; }),
                           found_roots_falsi.end());
    cout << "Metoda Regula Falsi (" << found_roots_falsi.size() << " unikalnych pierwiastkow):\n";
    for(double r : found_roots_falsi) cout << r << "\n";


    // Metoda Bisekcji
    for (double a_loop = start; a_loop < end; a_loop += STEP) {
        double b_loop = a_loop + STEP;
        if (b_loop > end) b_loop = end;

        double fa = f(a_loop), fb = f(b_loop);
        if (!isfinite(fa) || !isfinite(fb)) continue;

        if (fa * fb <= 0) {
            vector<double> iteracje;
            double root = metoda_bisekcji(f, a_loop, b_loop, iteracje);
            if (!isnan(root) && root >= a_loop - EPSILON && root <= b_loop + EPSILON) {
                found_roots_bisekcja.push_back(root);
            }
        }
        else if (std::fabs(fa) < EPSILON) {
             found_roots_bisekcja.push_back(a_loop);
        }
        if (a_loop + STEP >= end && std::fabs(fb) < EPSILON) {
             found_roots_bisekcja.push_back(b_loop);
        }
    }
    sort(found_roots_bisekcja.begin(), found_roots_bisekcja.end());
    found_roots_bisekcja.erase(unique(found_roots_bisekcja.begin(), found_roots_bisekcja.end(),
                                     [](double r1, double r2){ return std::fabs(r1-r2) < EPSILON*10; }),
                              found_roots_bisekcja.end());
    cout << "Metoda Bisekcji (" << found_roots_bisekcja.size() << " unikalnych pierwiastkow):\n";
    for(double r : found_roots_bisekcja) cout << r << "\n";


    // Metoda Siecznych
    for (double a_loop = start; a_loop < end; a_loop += STEP) {
        double b_loop = a_loop + STEP;
        if (b_loop > end) b_loop = end;

        double fa = f(a_loop), fb = f(b_loop);
        if (!isfinite(fa) || !isfinite(fb)) continue;

        if (fa * fb <= 0) {
            vector<double> iteracje;
            double root = metoda_siecznych(f, a_loop, b_loop, iteracje);
            if (!isnan(root) && root >= a_loop - EPSILON && root <= b_loop + EPSILON) {
                found_roots_sieczne.push_back(root);
            }
        }
        // Dodatkowe sprawdzenie punktów brzegowych, które mogłyby być pierwiastkami
        else if (std::fabs(fa) < EPSILON) {
             found_roots_sieczne.push_back(a_loop);
        }
        if (a_loop + STEP >= end && std::fabs(fb) < EPSILON) { // Dla ostatniego punktu 'b_loop'
             found_roots_sieczne.push_back(b_loop);
        }
    }
    std::ranges::sort(found_roots_sieczne.begin(), found_roots_sieczne.end());
    found_roots_sieczne.erase(unique(found_roots_sieczne.begin(), found_roots_sieczne.end(),
                                   [](double r1, double r2){ return std::fabs(r1-r2) < EPSILON*10; }),
                             found_roots_sieczne.end());
    cout << "Metoda Siecznych (" << found_roots_sieczne.size() << " unikalnych pierwiastkow):\n";
    for(double r : found_roots_sieczne) std::cout << r << "\n";
}
