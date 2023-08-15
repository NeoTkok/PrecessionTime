#ifndef TEST_CPP
#define TEST_CPP
 
#include <gtest/gtest.h>

#include "sitellite.hpp"
#include <Eigen/Dense>



constexpr double J2 = 1.082'626'68 * 1e-3;
constexpr double Re = 6'378'137;
constexpr double GM = 3.986'004'418'888'88 * 1e14;  



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
    for(double i = - M_PI; i < M_PI;i += M_PI / 15)
    {
        Keplerian A{10'000'000, 1./2., M_PI - 1, i};
        
        Cartesian B = toCartesian(A,GM);
        Keplerian C = toKeplerian(B,GM);
        
        ASSERT_NEAR(A.sem_axis- C.sem_axis,0., 1e-7);
        ASSERT_NEAR(A.eccentricity - C.eccentricity, 0., 1e-7);
        
        
        ASSERT_NEAR(A.anomaly - C.anomaly, 0. , 1e-10);
        ASSERT_NEAR(A.arg_peri - C.arg_peri, 0. , 1e-10);

    }
}

TEST(Cartesian, to_Cartesian_2) {
    for(double i = - M_PI; i < M_PI ;i += 0.1)
    {
        Keplerian A{10'000'000, 1./2., (M_PI - i >= M_PI)?(M_PI-i - 2*M_PI):(M_PI - i), i};
        Cartesian B = toCartesian(A,GM);
        Keplerian C = toKeplerian(B,GM);
        ASSERT_NEAR(A.sem_axis- C.sem_axis,0., 1e-3);
        ASSERT_NEAR(A.eccentricity - C.eccentricity, 0., 1e-11);
        ASSERT_NEAR(A.arg_peri, C.arg_peri, 1e-7);
        ASSERT_NEAR(A.anomaly - C.anomaly, 0., 1e-7);
    }
}


TEST(Cartesian, to_Keplerian) {
    for(double i = 1.0; i < 1000 ;i += 1.)

    {
        Eigen::Vector2d R{8'000'000., 0.};
        Eigen::Vector2d V{i, 7500.};
        Cartesian A{R,V};
        Keplerian B = toKeplerian(A,GM);
        Cartesian C = toCartesian(B,GM);
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
    Keplerian B = toKeplerian(A,GM);
    
    Eigen::Vector2d delta_V{0., 10.};
    Keplerian C =  maneuver(A,delta_V,GM);
    ASSERT_NEAR((C.sem_axis*(1+C.eccentricity)-B.sem_axis*(1+B.eccentricity))/1000, 63.7487116, 1e-8);
    
    delta_V[1] -= 20;

    Cartesian T = toCartesian(C,GM);
    C = maneuver(T,delta_V,GM);
    ASSERT_NEAR((B.sem_axis*(1+B.eccentricity)-C.sem_axis*(1+C.eccentricity))/1000, 0., 1e-1);
  
}


TEST(function, precession_rate) {
    double d = precession_rate(Re + 800'000, 0, 56. / 180. * M_PI, Re, J2, GM);
    ASSERT_NEAR(d, -7.44298 * 1e-7, 1e-12);

}

/*
TEST(Cartesian, T_maneuver) {
    Keplerian K{25'000'000., 0.8, 0, 3.7850937};


    //Coordinates R(7'500'000., 0.);
    //Coordinates V(2000, 5000.);
    //Cartesian A(R,V);
    Cartesian A = KeplerianToCartesian(K,GM);
    std::cout << "Vx = " <<  A.v.x << "    Vy = " << A.v.y << std::endl;
    double rate = 500;
    for (double j = 0.0; j < 2. * M_PI; j += 0.01){
        Coordinates delta_V(rate * cos(j), rate * sin(j));
        //std::cout << j  << "  " << T_maneuver(A, delta_V, 0.1, 50. * M_PI / 280.)/60./24. << std::endl;
    }    
}
*/


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}



#endif