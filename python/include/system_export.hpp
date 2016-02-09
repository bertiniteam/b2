//
//  system_export.hpp
//  Xcode_b2
//
//  Created by Collins, James B. on 2/9/16.
//  Copyright (c) 2016 West Texas A&M University. All rights reserved.
//

#ifndef Xcode_b2_system_export_hpp
#define Xcode_b2_system_export_hpp
#include <bertini2/system.hpp>



#include "python_common.hpp"


namespace bertini{
	namespace python{
		
		using namespace bertini;
		
		void ExportSystem();
		
		
		
		
		///////// NamedSymbol class(abstract) ////////////////
		template<typename NodeBaseT>
		class SystemVisitor: public def_visitor<SystemVisitor<NodeBaseT> >
		{
			friend class def_visitor_access;
			
		public:
			template<class PyClass>
			void visit(PyClass& cl) const;
			
		private:
			
		};
		
	}
}


#endif
