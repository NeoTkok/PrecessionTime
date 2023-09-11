#ifndef SATELLITE_HPP
#define SATELLITE_HPP

#include <Eigen/Dense>
#include <cmath>
#include <optional>
#include <limits>

using Vector2d = Eigen::Vector<double, 2>;

struct Earth{
    double J2 = 1.082'626'68 * 1e-3;
    double Re = 6'378'137;
    double GM = 3.986'004'418'888'88 * 1e14;  
};

struct Cartesian{
    Vector2d position;
    Vector2d velocity;
};

struct Keplerian{
    double sem_axis;
    double eccentricity;
    double arg_peri;
    double anomaly;
};


Cartesian toCartesian(const Keplerian& keplerian, const double gravParam) noexcept{
    const double param = keplerian.sem_axis * (1 - keplerian.eccentricity * keplerian.eccentricity);
    const double anomalyPlusw = keplerian.anomaly + keplerian.arg_peri;
    const double cos_anomaly = std::cos(keplerian.anomaly);
    const double sin_anomaly = std::sin(keplerian.anomaly);

    const double R = param / (1 + keplerian.eccentricity * cos_anomaly);
    const double x = R * cos_anomaly;
    const double y = R * sin_anomaly;

    const double muDivPsqrt = sqrt(gravParam / param);    
    const double Vx = -muDivPsqrt * sin_anomaly;
    const double Vy = muDivPsqrt * (keplerian.eccentricity + cos_anomaly);
    
    const double c = std::cos(keplerian.arg_peri);
    const double s = std::sin(keplerian.arg_peri);
    return {{x*c - y*s, y*c + x*s},{Vx*c - Vy*s, Vy*c + Vx*s}};
}


std::optional<Keplerian> toKeplerian(const Cartesian& cartesian, const double gravParam) noexcept{
    const double velSqr = cartesian.velocity.squaredNorm();
    const double vel = std::sqrt(velSqr);
    const double positionNorm = cartesian.position.norm();
    const double muDivR = gravParam / positionNorm;
    const Vector2d eccVector = ((velSqr - muDivR) * cartesian.position - cartesian.velocity.dot(cartesian.position) * cartesian.velocity) / gravParam;


    const double ecc = eccVector.norm();  // may be optional?
    if (ecc >= 1){
        return std::nullopt;
    }

    const double periapsisArgument = (ecc != 0) ? std::atan2(eccVector.y(), eccVector.x()) : 0;
    const Vector2d e1 = ecc != 0 ? eccVector / ecc : Vector2d{1, 0};
    const Vector2d e2(-e1.y(), e1.x());
    const double anomaly = std::atan2(cartesian.position.dot(e2), cartesian.position.dot(e1));
    const double semimajor = gravParam / (2 * muDivR - velSqr);

    return Keplerian{semimajor, ecc, periapsisArgument, anomaly};
}


Cartesian maneuver(const Cartesian& cartesian, const Eigen::Vector2d& delta_velocity){
    return {cartesian.position, cartesian.velocity + delta_velocity};
}

const double precession_rate(const double sem_axis, const double ecc, const double i,
                        const Earth& earth) noexcept{
    const double ReDivP = earth.Re / sem_axis / (1 - ecc * ecc);
    return -3 * ReDivP * ReDivP / sem_axis * std::sqrt(earth.GM/sem_axis) * earth.J2 * std::cos(i) / 2;
}

struct OrbitCartesian{
    Cartesian inPlane;
    double inclination;
};  

std::optional<double> T_maneuver(const Keplerian& keplerian, const OrbitCartesian& orbit, const Eigen::Vector2d& dv,
                        const double DeltaOmega, const Earth& earth) noexcept{
    
    const Cartesian newCartesian = maneuver(orbit.inPlane, dv);
    const std::optional<Keplerian> newKeplerian = toKeplerian(newCartesian, earth.GM);

    if (newKeplerian == std::nullopt){
        return std::nullopt;
    }

    const double OmegaDot = precession_rate(keplerian.sem_axis, keplerian.eccentricity, orbit.inclination, earth);
    const double newOmegaDot = precession_rate(newKeplerian->sem_axis, newKeplerian->eccentricity, orbit.inclination, earth);
    const double deltaOmegaDot = newOmegaDot - OmegaDot;
    
    const double k = (deltaOmegaDot * DeltaOmega < 0) ? deltaOmegaDot / std::abs(deltaOmegaDot) : 0;
    return (DeltaOmega + k * 2 * M_PI) / deltaOmegaDot;
}   

struct ReturnValue{
    double minTime = std::numeric_limits<double>::infinity();
    double anomaly = std::numeric_limits<double>::quiet_NaN();
    Vector2d dv = {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
};

struct AlgoParameters{
    int nAnomaly;
    int nAngle;
    int nVelocity;
    double maxVel;
};

struct OrbitKeplerian{
    Keplerian inPlane;
    double inclination;
};

std::vector<ReturnValue> time_precession(const OrbitKeplerian& orbitKepl, const double DeltaOmega,
                        const Earth& earth, const AlgoParameters& algParams) noexcept{
    const Cartesian cartesian = toCartesian(orbitKepl.inPlane, earth.GM); 
    std::vector<ReturnValue> result(algParams.nVelocity);

    const double delta_alpha = 2 * M_PI / algParams.nAnomaly;
    const double delta_anomaly = 2 * M_PI / algParams.nAngle;
    const double delta_vel = algParams.maxVel / algParams.nVelocity;

    for (unsigned int i_nu = 0; i_nu < algParams.nAnomaly; ++i_nu){ 
        const double nu = delta_anomaly * i_nu;
        const Keplerian newKeplerian = {orbitKepl.inPlane.sem_axis, orbitKepl.inPlane.eccentricity, orbitKepl.inPlane.arg_peri, nu};
        const Cartesian newCartesian = toCartesian(newKeplerian, earth.GM);
        for (unsigned int i_vel = 0; i_vel < algParams.nVelocity; ++i_vel){
            const double absV = delta_vel * i_vel;

            for (unsigned int i_alpha = 0; i_alpha < algParams.nAngle; ++i_alpha){
                const double alpha = i_alpha * delta_alpha;

                const Vector2d Vel = {absV * std::cos(alpha), absV * std::sin(alpha)};
                const std::optional<double> t = T_maneuver(newKeplerian, {newCartesian, orbitKepl.inclination}, Vel, DeltaOmega, earth); 
                
                if (t.has_value() && *t < result[i_vel].minTime){  
                    result[i_vel].minTime = *t; 
                    result[i_vel].anomaly = nu;
                    result[i_vel].dv = Vel;
                }

            }
        }    
    }
    return result;

}



#endif 
