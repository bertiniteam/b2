
 #include <iostream>
// after git changes remove lines 4 through 12 and add #include <bertini2/tracking.hpp> 
 #include <typeinfo>
 #include <boost/multiprecision/mpfr.hpp>
 #include <bertini2/limbo.hpp>
 #include <bertini2/mpfr_complex.hpp>
 #include <bertini2/start_system.hpp>
 #include <bertini2/tracking/tracker.hpp>
 #include <bertini2/num_traits.hpp>
 #include <bertini2/tracking/amp_cauchy_endgame.hpp>
 #include <bertini2/tracking/amp_powerseries_endgame.hpp>
 #include <bertini2/logging.hpp>

using namespace bertini::tracking;
using namespace bertini::tracking::endgame;

using TrackerType = AMPTracker; // select a tracker type
using CauchyEGType = EndgameSelector<TrackerType>::Cauchy;
using PSEGType = EndgameSelector<TrackerType>::PSEG;
using PrecisionConfig = TrackerTraits<TrackerType>::PrecisionConfig;
auto TestedPredictor = config::Predictor::HeunEuler;
unsigned ambient_precision = bertini::DoublePrecision();

using Variable = bertini::node::Variable;
using Var = std::shared_ptr<Variable>;
using VariableGroup = bertini::VariableGroup;

using System = bertini::System;

using dbl = std::complex<double>;
using mpfr = bertini::complex;
using mpfr_float = bertini::mpfr_float;
using mpq_rational = bertini::mpq_rational;

template<typename NumType> using Vec = Eigen::Matrix<NumType, Eigen::Dynamic, 1>;


using RealT = TrackerTraits<TrackerType>::BaseRealType; //change to RealT and ComplexT
using ComplexT = TrackerTraits<TrackerType>::BaseComplexType;

template<typename ...T>
ComplexT ComplexFromString(T... s)
{return bertini::NumTraits<ComplexT>::FromString(s...);}

template<typename ...T>
RealT RealFromString(T... s)
{return bertini::NumTraits<RealT>::FromString(s...);}

Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", " << ", ";");


bertini::LoggingInit(fatal);

    auto MakeVariable(std::string varname)
    {
        return std::make_shared<Variable>(varname);
    }

    auto Gamma()
    {
        return std::make_shared<bertini::node::Rational>(bertini::node::Rational::Rand());
    }

    void CauchyRun(bertini::System griewank_sys, CauchyEGType my_endgame, std::vector<Vec<ComplexT> > solutions_to_finish, ComplexT t_endgame_boundary);

    void PSEGRun(bertini::System griewank_sys, PSEGType my_endgame, std::vector<Vec<ComplexT> > solutions_to_finish, ComplexT t_endgame_boundary);


