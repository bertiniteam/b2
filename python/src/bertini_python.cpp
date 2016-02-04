
#include "bertini_python.hpp"


namespace bertini
{
	namespace python
	{
//		class Pure
//		{
//		public:
//			virtual void foo(unsigned int v)  = 0;
//			virtual void bar() = 0;
////			virtual int foo(double v) const = 0;
//		};
//		
//		class Empty : public Pure
//		{
//			
//		};
//		
//		
//		class Base : public virtual Empty
//		{
//		public:
//			
//			virtual void foo(unsigned int v) override
//			{
//				std::cout << v+5 << std::endl;
//			}
//				
//			
//		};
//				
//		
//		class Derived : public virtual Base
//		{
//		public:
//			Derived(float f)
//			{
//				num = f;
//			}
//			
//			
//			void bar() override
//			{
//				std::cout << "derived bar\n";
//			}
//			
//		private:
//			float num;
//		};
//		
//				
//				
//		struct PureWrap : Pure, wrapper<Pure>
//		{
//			
//			void foo(unsigned int v)
//			{
//				this->get_override("foo")(v);
//			}
//			
//		}; 
		


		BOOST_PYTHON_MODULE(libpybertini) // this name must match the name of the generated .so file.
		{
			
			/******** Testing(BEGIN) **************/
//			implicitly_convertible<std::shared_ptr<Derived>, std::shared_ptr<Pure> >();
			
//			class_<PureWrap,boost::noncopyable, std::shared_ptr<Pure> >("Pure", no_init)
//			.def("foo", pure_virtual(&Pure::foo) )
//			.def("bar", pure_virtual(&Pure::bar) )
//			;
//			
//			class_<Base, bases<Pure>, boost::noncopyable, std::shared_ptr<Base> >("Base", no_init )
//			;
//			
//			class_<Derived, bases<Base>, std::shared_ptr<Derived> >("Derived", init<float>() )
//			;
			
			
			
			/******** Testing(END) **************/
			
			
			ExportMpfr();
		
			SetupFunctionTree();
			ExportNode();
			ExportSymbols();
			//	ExportOperators();
			//	ExportRoots();
			//	ExportSystem();
			
		}
	
	}
}


