// python/function_tree.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with  python/function_tree.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
//  Jeb Collins
//  West Texas A&M University
//  Mathematics
//  Fall 2015
//
//
//  python/operator.cpp:  the source file for the python interface for operator classes.

#include <stdio.h>
#include "operator.hpp"

using namespace boost::python;
using namespace bertini::node;


namespace bertini{
	namespace python{
		
		

		
		
		void ExportOperators()
		{
			using Nodeptr = std::shared_ptr<node::Node>;
			
			class_< node::Operator, bases<node::Node>, boost::noncopyable >("Operator", no_init);
			
			// Unary Operators
			class_< node::UnaryOperator, bases<Operator>, boost::noncopyable >("UnaryOperator", no_init)
			.def("reset", &node::UnaryOperator::Reset)
			.def("set_child", &node::UnaryOperator::SetChild)
			.def("first_child", &node::UnaryOperator::first_child)
			.def("degree", unOpDeg1, unOpDeg1Overloads())
			.def("degree", unOpDeg2)
			.def("multidegree", &node::UnaryOperator::MultiDegree)
			.def("is_homogeneous", unOpIsHom1, unOpIsHom1Overloads())
			.def("is_homogeneous", unOpIsHom2)
			.def("homogenize", &node::UnaryOperator::Homogenize)
			.def("precision", &node::UnaryOperator::precision)
			;
			
			
			// Binary Operators
			class_< BinaryOperator, bases<Operator>, boost::noncopyable >("BinaryOperator", no_init)
			.def("reset", pure_virtual(&BinaryOperator::Reset))
			;
			
			
			// Nary Operators
			class_< node::NaryOperator, bases<Operator>, bases<node::Node>, boost::noncopyable >("NaryOperator", no_init)
			.def("reset", &node::NaryOperator::Reset)
			.def("add_child", &node::NaryOperator::AddChild)
			.def("first_child", &node::NaryOperator::first_child)
			.def("precision", &node::NaryOperator::precision)
			.def("children_size", &node::NaryOperator::children_size)
			;
			
			
			
			
			
			// Sum Operator
			class_<node::SumOperator, bases<node::NaryOperator>,std::shared_ptr<node::SumOperator> >
			("SumOperator", init<>())
			.def(init<const std::shared_ptr<Node> &,const std::shared_ptr<Node> &>())
			.def(init<const std::shared_ptr<Node> &,bool,const std::shared_ptr<Node> &,bool>())
			
			.def("add_child", sumAddChild1)
			.def("add_child", sumAddChild2)
			
			.def("degree", sumDeg1, sumDeg1Overloads())
			.def("degree", sumDeg2)
			.def("multidegree", &node::SumOperator::MultiDegree)
			.def("is_homogeneous", sumIsHom1, sumIsHom1Overloads())
			.def("is_homogeneous", sumIsHom2)
			.def("homogenize", &node::SumOperator::Homogenize)
			
			.def("differentiate", &node::SumOperator::Differentiate)
			
//			.def(self += other<Nodeptr>())
//			.def(self -= other<Nodeptr>())
			;
			
			
			// Negate Operator
			class_<node::NegateOperator, bases<node::UnaryOperator>,std::shared_ptr<node::NegateOperator> >
			("NegateOperator", init< optional<const Nodeptr&> >())
			.def("differentiate", &node::NegateOperator::Differentiate)
			.def("is_homogeneous", negIsHom1, negIsHom1Overloads())
			.def("is_homogeneous", negIsHom2)
			;
			
			
			// Multiplication Operator
			class_<node::MultOperator, bases<node::NaryOperator>,std::shared_ptr<node::MultOperator> >
			("MultOperator", init<>())
			.def(init<const std::shared_ptr<Node> &,const std::shared_ptr<Node> &>())
			.def(init<const std::shared_ptr<Node> &,bool,const std::shared_ptr<Node> &,bool>())
			
			.def("add_child", multAddChild1)
			.def("add_child", multAddChild2)
			
			.def("degree", multDeg1, multDeg1Overloads())
			.def("degree", multDeg2)
			.def("is_homogeneous", multIsHom1, multIsHom1Overloads())
			.def("is_homogeneous", multIsHom2)
			.def("homogenize", &node::MultOperator::Homogenize)
			
			.def("differentiate", &MultOperator::Differentiate)
			;
			
			
			// Power Operator(with any exponent)
			class_<node::PowerOperator, bases<BinaryOperator>, std::shared_ptr<node::PowerOperator> >
			("PowerOperator", init<>())
			.def(init<const std::shared_ptr<Node> &,const std::shared_ptr<Node> &>())
			
			.def("set_base", &node::PowerOperator::SetBase)
			.def("set_exponent", &node::PowerOperator::SetBase)
			
			.def("degree", powDeg1, powDeg1Overloads())
			.def("degree", powDeg2)
			.def("multidegree", &PowerOperator::MultiDegree)
			.def("is_homogeneous", powIsHom1, powIsHom1Overloads())
			.def("is_homogeneous", powIsHom2)
			.def("homogenize", &PowerOperator::Homogenize)
			.def("precision", &PowerOperator::precision)
			
			.def("reset", &PowerOperator::Reset)
			.def("differentiate", &PowerOperator::Differentiate)
			;
			
			
			// IntegerPower Operator(with integer exponents)
			class_<IntegerPowerOperator, bases<UnaryOperator>, std::shared_ptr<IntegerPowerOperator> >
			("IntegerPowerOperator", init<>())
			.def(init<const std::shared_ptr<Node> &,optional<int> >())
			
			.add_property("exponent", &IntegerPowerOperator::exponent, &IntegerPowerOperator::set_exponent)
			
			.def("degree", intpowDeg1, intpowDeg1Overloads())
			.def("is_homogeneous", intpowIsHom1, intpowIsHom1Overloads())
			.def("is_homogeneous", powIsHom2)
			;

			
			// Sqrt Operator
			class_<SqrtOperator, bases<UnaryOperator>, std::shared_ptr<SqrtOperator> >
			("SqrtOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("degree", sqrtDeg1, sqrtDeg1Overloads())
			
			.def("differentiate", &SqrtOperator::Differentiate)
			;

			
			// Exp Operator
			class_<ExpOperator, bases<UnaryOperator>, std::shared_ptr<ExpOperator> >
			("ExpOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("degree", expDeg1, expDeg1Overloads())
			
			.def("differentiate", &ExpOperator::Differentiate)
			;

			
			// Log Operator
			class_<LogOperator, bases<UnaryOperator>, std::shared_ptr<LogOperator> >
			("LogOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("degree", logDeg1, logDeg1Overloads())
			
			.def("differentiate", &LogOperator::Differentiate)
			;
			
			
			
			// TrigOperator
			class_<TrigOperator, bases<UnaryOperator>, boost::noncopyable >
			("TrigOperator", no_init)
			
			.def("degree", &TrigOperator::Degree)
			;

			
			
			// SinOperator
			class_<SinOperator, bases<TrigOperator>, std::shared_ptr<SinOperator> >
			("SinOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &SinOperator::Differentiate)
			;

			
			// CosOperator
			class_<CosOperator, bases<TrigOperator>, std::shared_ptr<CosOperator> >
			("CosOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &CosOperator::Differentiate)
			;

			
			// TanOperator
			class_<TanOperator, bases<TrigOperator>, std::shared_ptr<TanOperator> >
			("TanOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &TanOperator::Differentiate)
			;

			
			// ArcSinOperator
			class_<ArcSinOperator, bases<TrigOperator>, std::shared_ptr<ArcSinOperator> >
			("ArcSinOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &ArcSinOperator::Differentiate)
			;
			
			
			// ArcCosOperator
			class_<ArcCosOperator, bases<TrigOperator>, std::shared_ptr<ArcCosOperator> >
			("ArcCosOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &ArcCosOperator::Differentiate)
			;
			
			
			// ArcTanOperator
			class_<ArcTanOperator, bases<TrigOperator>, std::shared_ptr<ArcTanOperator> >
			("ArcTanOperator", init< optional<const std::shared_ptr<Node> &> >())
			
			.def("differentiate", &ArcTanOperator::Differentiate)
			;

		}
		
	} // re: python
} // re: bertini
