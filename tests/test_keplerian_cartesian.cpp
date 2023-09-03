#ifndef TEST_CPP
#define TEST_CPP
 
#include <gtest/gtest.h>

#include "sitellite.hpp"
#include <Eigen/Dense>


TEST(Keplerian, definition_1) {
    Keplerian A{10., 1./std::sqrt(2), 0.3, 3.};
    
    ASSERT_NEAR(A.sem_axis,10, 1e-15);
    ASSERT_NEAR(A.eccentricity, 1/sqrt(2), 1e-15);
    ASSERT_NEAR(A.arg_peri, 0.3, 1e-15);
    ASSERT_NEAR(A.anomaly, 3., 1e-15);
}

TEST(Keplerian, definition_2) {
    Keplerian A;

    ASSERT_NEAR(A.sem_axis,0., 1e-15);
    ASSERT_NEAR(A.eccentricity, 0., 1e-15);
    ASSERT_NEAR(A.arg_peri, 0., 1e-15);
    ASSERT_NEAR(A.anomaly,0., 1e-15);
}

TEST(Coordinates, definition) {
    Eigen::Vector2d A;
    Eigen::Vector2d B{23., -3};

    ASSERT_NEAR(A[0], 0, 1e-15);
    ASSERT_NEAR(A[1], 0, 1e-15);
    ASSERT_NEAR(B[0], 23., 1e-15);
    ASSERT_NEAR(B[1], -3., 1e-15);
}

TEST(Cartesian, definition_1) {
    Eigen::Vector2d R(10., 1./std::sqrt(2));
    Eigen::Vector2d V(0.3, 3.);

    Cartesian A{R,V};
    
    ASSERT_NEAR(A.position[0], 10, 1e-15);
    ASSERT_NEAR(A.position[1], 1/std::sqrt(2), 1e-15);
    ASSERT_NEAR(A.velocity[0], 0.3, 1e-15);
    ASSERT_NEAR(A.velocity[1], 3., 1e-15);
}

TEST(Cartesian, definition_2) {
    Cartesian A;
    
    ASSERT_NEAR(A.position[0], 0., 1e-15);
    ASSERT_NEAR(A.position[1], 0., 1e-15);
    ASSERT_NEAR(A.velocity[0], 0., 1e-15);
    ASSERT_NEAR(A.velocity[1], 0., 1e-15);
}

TEST(Cartesian, to_Cartesian_1) {
    Earth E;

    for(double i = - M_PI; i < M_PI;i += M_PI / 15)
    {
        Keplerian A{10'000'000, 1./2., M_PI - 1, i};
        
        Cartesian B = toCartesian(A,E.GM);
        std::optional<Keplerian> C = toKeplerian(B,E.GM);
        
        ASSERT_NEAR(A.sem_axis- C->sem_axis,0., 1e-7);
        ASSERT_NEAR(A.eccentricity - C->eccentricity, 0., 1e-7);
        
        
        ASSERT_NEAR(A.anomaly - C->anomaly, 0. , 1e-10);
        ASSERT_NEAR(A.arg_peri - C->arg_peri, 0. , 1e-10);

    }
}

TEST(Cartesian, to_Cartesian_2) {
    Earth E;
    for(double i = - M_PI; i < M_PI ;i += 0.1)
    {
        Keplerian A{10'000'000, 1./2., (M_PI - i >= M_PI)?(M_PI-i - 2*M_PI):(M_PI - i), i};
        Cartesian B = toCartesian(A,E.GM);
        std::optional<Keplerian> C = toKeplerian(B,E.GM);
        ASSERT_NEAR(A.sem_axis- C->sem_axis,0., 1e-3);
        ASSERT_NEAR(A.eccentricity - C->eccentricity, 0., 1e-11);
        ASSERT_NEAR(A.arg_peri, C->arg_peri, 1e-7);
        ASSERT_NEAR(A.anomaly - C->anomaly, 0., 1e-7);
    }
}


TEST(Cartesian, to_Keplerian) {
    Earth E;
    for(double i = 1.0; i < 1000 ;i += 1.){
        Eigen::Vector2d R{8'000'000., 0.};
        Eigen::Vector2d V{i, 7500.};
        Cartesian A{R,V};
        std::optional<Keplerian> B = toKeplerian(A,E.GM);
        Cartesian C = toCartesian(*B,E.GM);
        ASSERT_NEAR(A.position[0] - C.position[0], 0., 1e-8);
        ASSERT_NEAR(A.position[1] - C.position[1], 0., 1e-8);
        ASSERT_NEAR(A.velocity[0] - C.velocity[0], 0., 1e-5);
        ASSERT_NEAR(A.velocity[1] - C.velocity[1], 0., 1e-5);
    }
}

TEST(Cartesian, maneuver) {
    Eigen::Vector2d R{8'000'000., 0.};
    Eigen::Vector2d V{0., 7500.};
    Cartesian A{R,V};
    Earth E;

    std::optional<Keplerian> B = toKeplerian(A,E.GM);
    
    Eigen::Vector2d delta_V{0., 10.};
    std::optional<Keplerian> C = maneuver(A,delta_V,E.GM);
    ASSERT_NEAR((C->sem_axis*(1+C->eccentricity)-B->sem_axis*(1+B->eccentricity))/1000, 63.7487116, 1e-8);
    
    delta_V[1] -= 20;

    Cartesian T = toCartesian(*C,E.GM);
    C = maneuver(T,delta_V,E.GM);
    ASSERT_NEAR((B->sem_axis*(1+B->eccentricity)-C->sem_axis*(1+C->eccentricity))/1000, 0., 1e-1);
  
}

TEST(function, precession_rate) {
    Earth E;
    double d = precession_rate(E.Re + 800'000, 0, 56. / 180. * M_PI, E);
    ASSERT_NEAR(d, -7.44298 * 1e-7, 1e-12);

}



TEST(Cartesian, T_maneuver) {
    Keplerian K{6'800'000., 0., 0, 0};
    Earth E;

    Cartesian A = toCartesian(K,E.GM);

    const int nAng = 100;
    const int nAnom = 100;
    const int nVel = 100;
    const double Vmax = 100;

    std::vector<ReturnValue> X = time_precession({A, M_PI*97./180}, M_PI*5./180, E, {nAnom, nAng, nVel, Vmax}); 
    
   // for (int i = 0; i < nVel; ++i){
 //       std::cout << i * Vmax / nVel << "  " << X[i].minTime << std::endl;;
 //   }

}



int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif