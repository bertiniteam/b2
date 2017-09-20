// this file includes all cpp files for this project, and represents the Unity style build.

#include "src/basics/limbo.cpp"
#include "src/basics/mpfr_complex.cpp"
#include "src/basics/mpfr_extensions.cpp"

#include "src/function_tree/node.cpp"
#include "src/function_tree/operators/arithmetic.cpp"
#include "src/function_tree/operators/operator.cpp"
#include "src/function_tree/operators/trig.cpp"

#include "src/function_tree/roots/function.cpp"
#include "src/function_tree/roots/jacobian.cpp"

#include "src/function_tree/symbols/differential.cpp"
#include "src/function_tree/symbols/number.cpp"
#include "src/function_tree/symbols/special_number.cpp"
#include "src/function_tree/symbols/symbol.cpp"
#include "src/function_tree/symbols/variable.cpp"

#include "src/parallel/initialize_finalize.cpp"
#include "src/parallel/parallel.cpp"


#include "src/system/precon.cpp"
#include "src/system/slice.cpp"
#include "src/system/start_base.cpp"
#include "src/system/system.cpp"

#include "src/system/start/total_degree.cpp"
#include "src/system/start/mhom.cpp"
#include "src/system/start/user.cpp"

#include "src/tracking/explicit_predictors.cpp"

