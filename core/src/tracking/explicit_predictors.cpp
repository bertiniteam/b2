//This file is part of Bertini 2.0.
//
//heun_euler.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//heun_euler.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with heun_euler.hpp.  If not, see <http://www.gnu.org/licenses/>.
//

//  predictor.hpp
//
//  copyright 2015
//  James B. Collins
//  West Texas A&M University
//  Department of Mathematics
//  Spring 2016


/**
 \file explicit_predictors.hpp
 
 \brief Defines all static variables for ExplicitRKPredictor class
 */

#include "bertini2/trackers/explicit_predictors.hpp"

namespace bertini{
	namespace tracking{
		namespace predict{


// Ptr must be filled row first.
// i.e. aPtr[] = {a_11, a_21, a_12, a_22}

/* Euler Butcher Table
 |
 0|0
 -----
 |1
 */


const mpq_rational ExplicitRKPredictor::aEulerPtr_[] = {mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor::aEuler_(aEulerPtr_);

const mpq_rational ExplicitRKPredictor::bEulerPtr_[] = {mpq_rational(1,1)};

const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor::bEuler_(bEulerPtr_);

const mpq_rational ExplicitRKPredictor::cEulerPtr_[] = {mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,1,1> ExplicitRKPredictor::cEuler_(cEulerPtr_);



/* Heun-Euler Butcher Table
 0 |
 1 | 1
 -----------
 |1/2  1/2
 | 1    0
 */


const mpq_rational ExplicitRKPredictor::aHeunEulerPtr_[] = {mpq_rational(0,1), mpq_rational(1,1),
	mpq_rational(0,1), mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,2,2> ExplicitRKPredictor::aHeunEuler_(aHeunEulerPtr_);

const mpq_rational ExplicitRKPredictor::bHeunEulerPtr_[] = {mpq_rational(1,2), mpq_rational(1,2)};

const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor::bHeunEuler_(bHeunEulerPtr_);

const mpq_rational ExplicitRKPredictor::b_minus_bstarHeunEulerPtr_[] = {mpq_rational(-1,2), mpq_rational(1,2)};

const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor::b_minus_bstarHeunEuler_(b_minus_bstarHeunEulerPtr_);

const mpq_rational ExplicitRKPredictor::cHeunEulerPtr_[] = {mpq_rational(0,1), mpq_rational(1,1)};

const Eigen::Matrix<mpq_rational,2,1> ExplicitRKPredictor::cHeunEuler_(cHeunEulerPtr_);



/* RK4 Butcher Table
 0   | 0   0  0  0
 1/2 | 1/2 0  0  0
 1/2 | 0  1/2 0  0
 1   | 0   0  1  0
 ------------------
 | 1/6 1/3 1/3 1/6
 */


const mpq_rational ExplicitRKPredictor::aRK4Ptr_[] = {mpq_rational(0,1), mpq_rational(1,2),
	mpq_rational(0,1), mpq_rational(0,1),
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1,2), mpq_rational(0,1),
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1,1),
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,4,4> ExplicitRKPredictor::aRK4_(aRK4Ptr_);

const mpq_rational ExplicitRKPredictor::bRK4Ptr_[] = {mpq_rational(1,6), mpq_rational(1,3),
	mpq_rational(1,3),mpq_rational(1,6)};

const Eigen::Matrix<mpq_rational,4,1> ExplicitRKPredictor::bRK4_(bRK4Ptr_);

const mpq_rational ExplicitRKPredictor::cRK4Ptr_[] = {mpq_rational(0,1), mpq_rational(1,2),
	mpq_rational(1,2),mpq_rational(1,1)};

const Eigen::Matrix<mpq_rational,4,1> ExplicitRKPredictor::cRK4_(cRK4Ptr_);




/* RKF45 Butcher Table
 https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 */


const mpq_rational ExplicitRKPredictor::aRKF45Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,4),mpq_rational(3,32), mpq_rational(1932,2197),mpq_rational(439,216),mpq_rational(-8,27),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(9,32), mpq_rational(-7200,2197),mpq_rational(-8,1),mpq_rational(2,1),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(7296,2197),mpq_rational(3680,513),mpq_rational(-3544,2565),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(-845,4104),mpq_rational(1859,4104),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(-11,40),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,6,6> ExplicitRKPredictor::aRKF45_(aRKF45Ptr_);


const mpq_rational ExplicitRKPredictor::bRKF45Ptr_[] =
{mpq_rational(16,135), mpq_rational(0,1), mpq_rational(6656,12825), mpq_rational(28561,56430),mpq_rational(-9,50),mpq_rational(2,55)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::bRKF45_(bRKF45Ptr_);


const mpq_rational ExplicitRKPredictor::b_minus_bstarRKF45Ptr_[] =
{mpq_rational(1,360), mpq_rational(0,1), mpq_rational(-128,4275), mpq_rational(-2197,75240),mpq_rational(1,50),mpq_rational(2,55)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::b_minus_bstarRKF45_(b_minus_bstarRKF45Ptr_);


const mpq_rational ExplicitRKPredictor::cRKF45Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,4), mpq_rational(3,8), mpq_rational(12,13),mpq_rational(1,1),mpq_rational(1,2)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::cRKF45_(cRKF45Ptr_);




/* RK Cash-Karp45 Butcher Table
 https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
 */


const mpq_rational ExplicitRKPredictor::aRKCK45Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,5),mpq_rational(3,40), mpq_rational(3,10),mpq_rational(-11,54),mpq_rational(1631,55296),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(9,40), mpq_rational(-9,10),mpq_rational(5,2),mpq_rational(175,512),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(6,5),mpq_rational(-70,27),mpq_rational(575,13824),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(35,27),mpq_rational(44275,110592),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(253,4096),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,6,6> ExplicitRKPredictor::aRKCK45_(aRKCK45Ptr_);


const mpq_rational ExplicitRKPredictor::bRKCK45Ptr_[] =
{mpq_rational(37,378), mpq_rational(0,1), mpq_rational(250,621), mpq_rational(125,594),mpq_rational(0,1),mpq_rational(512,1771)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::bRKCK45_(bRKCK45Ptr_);


const mpq_rational ExplicitRKPredictor::b_minus_bstarRKCK45Ptr_[] =
{mpq_rational(-277,64512), mpq_rational(0,1), mpq_rational(6925,370944), mpq_rational(-6925,202752),mpq_rational(-277,14336),mpq_rational(277,7084)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::b_minus_bstarRKCK45_(b_minus_bstarRKCK45Ptr_);


const mpq_rational ExplicitRKPredictor::cRKCK45Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,5), mpq_rational(3,10), mpq_rational(3,5),mpq_rational(1,1),mpq_rational(7,8)};

const Eigen::Matrix<mpq_rational,6,1> ExplicitRKPredictor::cRKCK45_(cRKCK45Ptr_);




/* RK Dormand-Prince56 Butcher Table
 Prince and Dormand.  High order embedded Runge-Kutta formulae.  J. Comput. Appl. Math., 7(1):67-75, 1981
 */


const mpq_rational ExplicitRKPredictor::aRKDP56Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,10),mpq_rational(-2,81), mpq_rational(615,1372),mpq_rational(3243,5500),mpq_rational(-26492,37125),mpq_rational(5561,2376),mpq_rational(465467,266112),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(20,81), mpq_rational(-270,343),mpq_rational(-54,55),mpq_rational(72,55),mpq_rational(-35,11),mpq_rational(-2945,1232),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1053,1372),mpq_rational(50949,71500),mpq_rational(2808,23375),mpq_rational(-24117,31603),mpq_rational(-5610201,14158144),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(4998,17875),mpq_rational(-24206,37125),mpq_rational(899983,200772),mpq_rational(10513573,3212352),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(338,459),mpq_rational(-5225,1836),mpq_rational(-424325,205632),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(3925,4056),mpq_rational(376225,454272),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,8,8> ExplicitRKPredictor::aRKDP56_(aRKDP56Ptr_);


const mpq_rational ExplicitRKPredictor::bRKDP56Ptr_[] =
{mpq_rational(821,10800), mpq_rational(0,1), mpq_rational(19683,71825), mpq_rational(175273,912600),mpq_rational(395,3672),mpq_rational(785,2704),mpq_rational(3,50),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor::bRKDP56_(bRKDP56Ptr_);


const mpq_rational ExplicitRKPredictor::b_minus_bstarRKDP56Ptr_[] =
{mpq_rational(13,2400), mpq_rational(0,1), mpq_rational(-19683,618800), mpq_rational(2401,31200),mpq_rational(-65,816),mpq_rational(15,416),mpq_rational(521,5600),mpq_rational(-1,10)};

const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor::b_minus_bstarRKDP56_(b_minus_bstarRKDP56Ptr_);


const mpq_rational ExplicitRKPredictor::cRKDP56Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,10), mpq_rational(2,9), mpq_rational(3,7),mpq_rational(3,5),mpq_rational(4,5),mpq_rational(1,1),mpq_rational(1,1)};

const Eigen::Matrix<mpq_rational,8,1> ExplicitRKPredictor::cRKDP56_(cRKDP56Ptr_);





/* RK Verner 67 Butcher Table
 J.T. Verner.  Explicit Runge-Kutta Methods with Estimates of the Local Truncation Error.  SIAM J. Num. Anal., 15(4):pp.772-790, 1978
 */


const mpq_rational ExplicitRKPredictor::aRKV67Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,12),mpq_rational(0,1), mpq_rational(1,16),mpq_rational(21,16),mpq_rational(1344688,250563),mpq_rational(-559,384),mpq_rational(-625,224),mpq_rational(-12253,99144),mpq_rational(30517,2512),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1,6), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(3,16),mpq_rational(-81,16),mpq_rational(-1709184,83521),mpq_rational(6,1),mpq_rational(12,1),mpq_rational(16,27),mpq_rational(-7296,157),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(9,2),mpq_rational(1365632,83521),mpq_rational(-204,47),mpq_rational(-456,47),mpq_rational(16,459),mpq_rational(268728,7379),
	
	mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(0,1),mpq_rational(0,1),mpq_rational(-78208,250563),mpq_rational(14,39),mpq_rational(48,91),mpq_rational(29072,161109),mpq_rational(2472,2041),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(-4913,78208),mpq_rational(14739,136864),mpq_rational(-2023,75816),mpq_rational(-3522621,10743824),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(6,7),mpq_rational(112,12393),mpq_rational(132,157),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(-12393,4396),
	
	mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1),mpq_rational(0,1)};

const Eigen::Matrix<mpq_rational,10,10> ExplicitRKPredictor::aRKV67_(aRKV67Ptr_);


const mpq_rational ExplicitRKPredictor::bRKV67Ptr_[] =
{mpq_rational(2881,40320), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(1216,2961),mpq_rational(-2624,4095),mpq_rational(24137569,57482880),mpq_rational(-4,21),mpq_rational(0,1),mpq_rational(4131,3920),mpq_rational(-157,1260)};

const Eigen::Matrix<mpq_rational,10,1> ExplicitRKPredictor::bRKV67_(bRKV67Ptr_);


const mpq_rational ExplicitRKPredictor::b_minus_bstarRKV67Ptr_[] =
{mpq_rational(-17,2688), mpq_rational(0,1), mpq_rational(0,1), mpq_rational(272,4935),mpq_rational(-272,273),mpq_rational(24137569,57482880),mpq_rational(-34,105),mpq_rational(-7,90),mpq_rational(4131,3920),mpq_rational(-157,1260)};

const Eigen::Matrix<mpq_rational,10,1> ExplicitRKPredictor::b_minus_bstarRKV67_(b_minus_bstarRKV67Ptr_);


const mpq_rational ExplicitRKPredictor::cRKV67Ptr_[] =
{mpq_rational(0,1), mpq_rational(1,12), mpq_rational(1,6), mpq_rational(1,4),mpq_rational(3,4),mpq_rational(16,17),mpq_rational(1,2),mpq_rational(1,1),mpq_rational(2,3),mpq_rational(1,1)};

const Eigen::Matrix<mpq_rational,10,1> ExplicitRKPredictor::cRKV67_(cRKV67Ptr_);



			
		} // re: predict
	}// re: tracking
}// re: bertini





