#ifndef SITELLITE_HPP
#define SITELLITE_HPP

#include <cmath>

constexpr double J2 = 1.082'626'68 * 1e-3;
constexpr double Re = 6'378'137;
constexpr double GM = 3.986'004'418'888'88 * 1e14;  


struct Coordinates{
    double x;
    double y;

    Coordinates() : x(0), y(0) {}

    Coordinates(const double x, const double y) : x(x), y(y){}
};


struct Cartesian{
    Coordinates r;
    Coordinates v;

    Cartesian(const Coordinates r, const Coordinates v) : r(r), v(v) {}
    
    Cartesian() : r(), v() {}
   };

// Перегрузка оператора << для структуры Coordinates
std::ostream& operator<<(std::ostream& os, const Cartesian& C) {
    os << "Rx = " << C.r.x << "\nRy = " << C.r.y << "\nVx = " << C.v.x << "\nRy = " << C.v.y << std::endl;
    return os;
}

struct Keplerian{
    double a;
    double e;
    double w;
    double nu;

    Keplerian(const double a, const double e, const double w,
            const double nu) : a(a), e(e), w(w), nu(nu) {}
    
    Keplerian() : a(0), e(0), w(0), nu(0) {}
};

std::ostream& operator<<(std::ostream& os, const Keplerian& K) {
    os << "a = " << K.a << "\ne = " << K.e << "\nw = " << K.w << "\nnu = " << K.nu << std::endl;
    return os;
}

Cartesian KeplerianToCartesian(Keplerian const K, double const g){
    Cartesian C;
    double phi = K.w + K.nu;

    C.r.x = K.a * (1 - K.e * K.e) * cos(phi) / (1. + K.e * cos(K.nu));
    C.r.y = K.a * (1 - K.e * K.e) * sin(phi) / (1. + K.e * cos(K.nu));

    double r = sqrt(C.r.x * C.r.x + C.r.y * C.r.y);
    double v = sqrt(g * (2. / r - 1. / K.a));
    double h = sqrt(g * K.a * (1 - K.e * K.e));
    double fpa = acos(h / (r * v));
    if (K.nu > M_PI)
        fpa = 2 * M_PI - fpa;
    
    C.v.x = v * sin(fpa - phi);
    C.v.y = v * cos(phi - fpa);

    return C;
}

Keplerian CartesianToKeplerian(Cartesian const C, double const g){
    Keplerian K;
    double h = C.r.x * C.v.y - C.r.y * C.v.x;
    double r = sqrt(C.r.x * C.r.x + C.r.y * C.r.y);
    double v = sqrt(C.v.x * C.v.x + C.v.y * C.v.y);
    K.a = 1. / (2. / r - v * v / g);
    K.e = sqrt(1 - h * h / g / K.a);
    double alpha = atan2(C.r.y, C.r.x) + 2 * M_PI; // = w + nu

    K.nu = acos((K.a * (1 - K.e * K.e) / r - 1) / K.e);
    
    K.w = alpha - K.nu;
    
    Cartesian new_C = KeplerianToCartesian(K,g); 
    double deltaV = (C.v.x - new_C.v.x)*(C.v.x - new_C.v.x) + (C.v.y - new_C.v.y)*(C.v.y - new_C.v.y);
    double deltaR = (C.r.x - new_C.r.x)*(C.r.x - new_C.r.x) + (C.r.y - new_C.r.y)*(C.r.y - new_C.r.y);
    
    if (deltaV + deltaR + 1 > 1) // чтоб убрать машинную неточность и повернуть истиную аномалию за афелий
        K.nu = 2 * M_PI - K.nu;
        
    K.w = (alpha <= K.nu) ? alpha - K.nu + 2 * M_PI: alpha - K.nu;

    return K;
}


Keplerian maneuver(const Cartesian C, const Coordinates dv, const double g){
    Cartesian new_C = C;
    new_C.v.x += dv.x;
    new_C.v.y += dv.y;
    
    Keplerian K = CartesianToKeplerian(new_C, g);
    
    return K;
}

double precession_rate(const double a, const double e, const double i) {
    return -3. * Re * Re / (a * a * a * sqrt(a) * (1 - e * e) * (1 - e * e)) *sqrt(GM)* J2 * cos(i) / 2.;
}

double T_maneuver(const Cartesian C, const Coordinates dv, const double DeltaOmega, const double i){
    Keplerian K = CartesianToKeplerian(C,GM);
    double omegaDot = precession_rate(K.a, K.e, i);
    
    Keplerian new_K = maneuver(C, dv, GM);
    double new_omegaDot = precession_rate(new_K.a, new_K.e, i);
    

    double delta_omegaDot = new_omegaDot - omegaDot;
    if (delta_omegaDot * DeltaOmega > 0)
        return DeltaOmega / delta_omegaDot;
    
    if (delta_omegaDot > 0 && DeltaOmega < 0)
        return (2 * M_PI + DeltaOmega) / delta_omegaDot;
    if (delta_omegaDot < 0 && DeltaOmega > 0)
        return (- 2 * M_PI + DeltaOmega) / delta_omegaDot;

    return 0.;
}


#endif 