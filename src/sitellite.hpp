#ifndef SITELLITE_HPP
#define SITELLITE_HPP

#include <Eigen/Dense>
#include <cmath>

constexpr double J2 = 1.082'626'68 * 1e-3;
constexpr double Re = 6'378'137;
constexpr double GM = 3.986'004'418'888'88 * 1e14;  


struct Cartesian{
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
};


struct Keplerian{
    double a;
    double e;
    double w;
    double nu;
};


Cartesian KeplerianToCartesian(const Keplerian& K, const double g){
    double phi = K.w + K.nu;
    Eigen::Vector3d r;
    Eigen::Vector3d v;

    r[0] = K.a * (1 - K.e * K.e) * std::cos(phi) / (1. + K.e * std::cos(K.nu));
    r[1] = K.a * (1 - K.e * K.e) * std::sin(phi) / (1. + K.e * std::cos(K.nu));

    double R = r.norm();
    double V = std::sqrt(g * (2. / R - 1. / K.a));
    double h = std::sqrt(g * K.a * (1 - K.e * K.e));
    double fpa = acos(h / (R * V));
    if (K.nu > M_PI)
        fpa = 2 * M_PI - fpa;
    
    v[0] = V * sin(fpa - phi);
    v[1] = V * cos(phi - fpa);

    return Cartesian{r,v};
}

Keplerian CartesianToKeplerian(const Cartesian& C, const double g){

    double h = C.position[0] * C.velocity[1] - C.position[1] * C.velocity[0];
    double R = C.position.norm();
    double V = C.velocity.norm();

    double a = 1. / (2. / R - V * V / g);
    double e = std::sqrt(1. - h * h / g / a);
    double alpha = atan2(C.position[1], C.position[0]) + 2 * M_PI; // = w + nu
    
    double nu = acos((a * (1 - e * e) / R - 1) / e);
    double w = alpha - nu;
    Keplerian K{a,e,w,nu};
    Cartesian new_C = KeplerianToCartesian(K,g); 
    double delta =  (C.velocity - new_C.velocity).norm() + (C.position - new_C.position).norm();
    
    if (delta > 1e-7)  
        K.nu = 2 * M_PI - K.nu;
        
    K.w = (alpha <= K.nu) ? alpha - K.nu + 2 * M_PI: alpha - K.nu;

    return K;
}


const Keplerian maneuver(const Cartesian& C, const Eigen::Vector3d& dv, const double g){
    Cartesian new_C = C;
    new_C.velocity += dv;
    
    Keplerian K = CartesianToKeplerian(new_C, g);
    
    return K;
}

const double precession_rate(const double a, const double e, const double i) {
    return -3. * Re*Re / (a*a*a * std::sqrt(a) * (1-e*e)*(1-e*e)) * std::sqrt(GM) * J2 * std::cos(i) / 2.;
}

double T_maneuver(const Cartesian& C, const Eigen::Vector3d& dv, const double DeltaOmega, const double i){
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