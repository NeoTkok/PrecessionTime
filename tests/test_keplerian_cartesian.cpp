#ifndef TEST_CPP
#define TEST_CPP

#include <gtest/gtest.h>

#include "sitellite.hpp"
#include <iostream>
#include <Eigen/Dense>


TEST(Keplerian, definition_1) {
    Keplerian A{10., 1./std::sqrt(2), 0.3, 3.};
    
    ASSERT_NEAR(A.a, 10, 1e-15);
    ASSERT_NEAR(A.e, 1/sqrt(2), 1e-15);
    ASSERT_NEAR(A.w, 0.3, 1e-15);
    ASSERT_NEAR(A.nu, 3., 1e-15);
}

TEST(Keplerian, definition_2) {
    Keplerian A;

    ASSERT_NEAR(A.a, 0., 1e-15);
    ASSERT_NEAR(A.e, 0., 1e-15);
    ASSERT_NEAR(A.w, 0., 1e-15);
    ASSERT_NEAR(A.nu,0., 1e-15);
}

TEST(Coordinates, definition) {
    Eigen::Vector3d A;
    Eigen::Vector3d B{23., -3, 0.};

    ASSERT_NEAR(A[0], 0, 1e-15);
    ASSERT_NEAR(A[1], 0, 1e-15);
    ASSERT_NEAR(B[0], 23., 1e-15);
    ASSERT_NEAR(B[1], -3., 1e-15);
}

TEST(Cartesian, definition_1) {
    Eigen::Vector3d R(10., 1./std::sqrt(2), 0.);
    Eigen::Vector3d V(0.3, 3., 0.);

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
    for(double i = 0.; i < M_PI*2 ;i += 0.1)
    {
        Keplerian A{10'000'000, 1./2., M_PI + 1, i};
        Cartesian B = KeplerianToCartesian(A, GM);
        Keplerian C = CartesianToKeplerian(B,GM);
        std::cout << "";

        ASSERT_NEAR(A.a - C.a, 0., 1e-7);
        ASSERT_NEAR(A.e - C.e, 0., 1e-7);
        ASSERT_NEAR(A.w - C.w, 0., 1e-1);
        ASSERT_NEAR(A.nu - C.nu, 0., 1e-5);
    }
}

TEST(Cartesian, to_Cartesian_2) {
    for(double i = 0.5; i < M_PI*2 ;i += 0.1)
    {
        Keplerian A{10'000'000, 1./2., 2 * M_PI - i, i};
        Cartesian B = KeplerianToCartesian(A, GM);
        Keplerian C = CartesianToKeplerian(B,GM);
        std::cout <<  "";

        ASSERT_NEAR(A.a - C.a, 0., 1e-3);
        ASSERT_NEAR(A.e - C.e, 0., 1e-11);
        ASSERT_NEAR(A.w, C.w, 1e-7);
        ASSERT_NEAR(A.nu - C.nu, 0., 1e-7);
    }
}


TEST(Cartesian, to_Keplerian) {
    for(double i = 1.0; i < 1000 ;i += 1.)

    {
        Eigen::Vector3d R{8'000'000., 0., 0.};
        Eigen::Vector3d V{i, 7500., 0.};
        Cartesian A{R,V};
        Keplerian B = CartesianToKeplerian(A,GM);
        Cartesian C = KeplerianToCartesian(B,GM);
        ASSERT_NEAR(A.position[0] - C.position[0], 0., 1e-8);
        ASSERT_NEAR(A.position[1] - C.position[1], 0., 1e-8);
        ASSERT_NEAR(A.velocity[0] - C.velocity[0], 0., 1e-5);
        ASSERT_NEAR(A.velocity[1] - C.velocity[1], 0., 1e-5);
    }
}

TEST(Cartesian, maneuver) {
    Eigen::Vector3d R{8'000'000., 0., 0.};
    Eigen::Vector3d V{0., 7500., 0.};
    Cartesian A{R,V};
    Keplerian B = CartesianToKeplerian(A,GM);
    
    Eigen::Vector3d delta_V{0., 10. ,0.};
    Keplerian C =  maneuver(A,delta_V,GM);
    ASSERT_NEAR((C.a*(1+C.e)-B.a*(1+B.e))/1000 , 63.7487116, 1e-8);
    
    delta_V[1] -= 20;
    std::cout << "";

    Cartesian T = KeplerianToCartesian(C,GM);
    C = maneuver(T,delta_V,GM);
    ASSERT_NEAR((B.a*(1+B.e)-C.a*(1+C.e))/1000, 0., 1e-1);
  
}


TEST(function, precession_rate) {
    double d = precession_rate(Re + 800'000, 0, 56. / 180. * M_PI);
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