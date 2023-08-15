#ifndef SITELLITE_HPP
#define SITELLITE_HPP

#include <Eigen/Dense>
#include <cmath>

using Vector2d = Eigen::Vector<double, 2>;

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


Keplerian toKeplerian(const Cartesian& cartesian, const double gravParam) noexcept {
    const double velSqr = cartesian.velocity.squaredNorm();
	const double vel = std::sqrt(velSqr);
	const double positionNorm = cartesian.position.norm();
	const double muDivR = gravParam / positionNorm;
	const Vector2d eccVector = ((velSqr - muDivR) * cartesian.position - cartesian.velocity.dot(cartesian.position) * cartesian.velocity) / gravParam;
	const double ecc = eccVector.norm();
	const double periapsisArgument = (ecc != 0) ? std::atan2(eccVector.y(), eccVector.x()) : 0;

	const Vector2d e1 = ecc != 0 ? eccVector / ecc : Vector2d{1, 0};
	const Vector2d e2(-e1.y(), e1.x());
	const double anomalyly = std::atan2(cartesian.position.dot(e2), cartesian.position.dot(e1));
	const double semimajor = gravParam / (2 * muDivR - velSqr);
    
	return {semimajor, ecc, periapsisArgument, anomalyly};
}


const Keplerian maneuver(const Cartesian& cartesian, const Eigen::Vector2d& delta_velocity, const double gravParam) noexcept{
    Cartesian new_cartesian = cartesian;
    new_cartesian.velocity += delta_velocity;    
    
    return toKeplerian(new_cartesian, gravParam);
}


const double precession_rate(const double sem_axis, const double ecc, const double i,
                        const double Re, const double J2, const double gravParam) noexcept{
    const double ReDivP = Re / sem_axis / (1 - ecc * ecc);
    return -3. * ReDivP * ReDivP / sem_axis * std::sqrt(gravParam/sem_axis) * J2 * std::cos(i) / 2.;
}


double T_maneuver(const Cartesian& cartesian, const Eigen::Vector2d& dv, const double DeltaOmega,
                const double i, const double Re, const double J2, const double gravParam) noexcept{
    const Keplerian keplerian = toKeplerian(cartesian,gravParam);
    const double omegaDot = precession_rate(keplerian.sem_axis, keplerian.eccentricity, i, Re, J2, gravParam);
    
    const Keplerian new_keplerian = maneuver(cartesian, dv, gravParam);
    const double new_omegaDot = precession_rate(new_keplerian.sem_axis, new_keplerian.eccentricity, i, Re, J2, gravParam);
        
    const double delta_omegaDot = new_omegaDot - omegaDot;
    const double k = (delta_omegaDot * DeltaOmega < 0) ? delta_omegaDot / std::abs(delta_omegaDot) : 0;

    return (delta_omegaDot != 0) ? (DeltaOmega + k * 2 * M_PI) / delta_omegaDot : -1;
}


#endif 