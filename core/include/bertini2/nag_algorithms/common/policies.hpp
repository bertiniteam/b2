//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/zero_dim_solve/policies.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/zero_dim_solve/policies.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/zero_dim_solve/policies.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) 2015 - 2017 by Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.


/**
\file bertini2/nag_algorithms/zero_dim_solve/policies.hpp 

\brief Provides some policies for the zero dim algorithm.

You can also provide your own, that's the point of these policies.
*/

#pragma once


namespace bertini {

	namespace policy{



		/**
		\brief A base class for system management for the zero-dim algorithm.
		*/
		template<	 typename D,
					 typename SystemType, typename StartSystemType
					,typename StoredSystemType, typename StoredStartSystemType>
		struct SysMgmtPolicy
		{

			using SystemT = SystemType;
			using StartSystemT = StartSystemType;

			using StoredSystemT = StoredSystemType;
			using StoredStartSystemT = StoredStartSystemType;

private:
			// convert the base endgame into the derived type.
			const D& AsDerived() const
			{
				return static_cast<const D&>(*this);
			}

			// convert the base endgame into the derived type.
			D& AsDerived() 			
			{
				return static_cast<D&>(*this);
			}

public:
			/**
			A getter for the system to be tracked to.
			*/
			const SystemT& TargetSystem() const
			{
				return AsDerived().target_system_;
			}

			/**
			A getter for the homotopy being used.
			*/
			const SystemT& Homotopy() const
			{
				return AsDerived().homotopy_;
			}

			/**
			A getter for the start system being used.
			*/
			const StartSystemT & StartSystem() const
			{
				return AsDerived().start_system_;
			}



			/**
			A setter for the system to be tracked to.
			*/
			void TargetSystem(StoredSystemT const& sys)
			{
				AsDerived().target_system_ = sys;
			}

			/**
			A setter for the homotopy being used.
			*/
			void Homotopy(StoredSystemT const& sys)
			{
				AsDerived().homotopy_ = sys;
			}

			/**
			A setter for the start system being used.
			*/
			void StartSystem(StoredStartSystemT const& sys)
			{
				AsDerived().start_system_ = sys;
			}
		};




		/**
		This system management policy makes it so that the zero dim algorithm makes a clone of the supplied system when the algorithm is created.  The zerodim algorithm, and others, will homogenize the system (that's why you want a clone, so it leaves your original system untouched), form a start system (of your type, inferred from the template parameter for the zerodim alg), and couple the two together into the homotopy used to track.

		If you don't want it to take copies, or homogenize, etc, use a different policy.

		\see RefToGiven
		*/
		template<typename SystemType, typename StartSystemType>
		struct CloneGiven : public SysMgmtPolicy<CloneGiven<SystemType, StartSystemType>, SystemType, StartSystemType, SystemType, StartSystemType>
		{

			using SMP = SysMgmtPolicy<CloneGiven<SystemType, StartSystemType>, SystemType, StartSystemType, SystemType, StartSystemType>;
			friend SMP;

			using StoredSystemT = typename SMP::StoredSystemT;
			using StoredStartSystemT = typename SMP::StoredStartSystemT;

			using SMP::TargetSystem;
			using SMP::StartSystem;
			using SMP::Homotopy;

			using SystemT = SystemType;
			using StartSystemT = StartSystemType;
			
private:
			StoredSystemT target_system_;
			StoredStartSystemT start_system_;
			StoredSystemT homotopy_;


public:
			/**
			Simply forward on the systems for the constructor.  The AtConstruct function is to be called by the user of this policy, at construct time.
			*/
			CloneGiven(SystemType const& target) : target_system_(AtConstruct(target))
			{}


			static
			StoredSystemT AtConstruct(SystemType const& sys)
			{
				return Clone(sys);
			}


			/**
			In contrast at AtConstruct, the AtSet function copies the given system into the stored system when setting, after construction.
			*/
			template<typename T>
			static
			T AtSet(T const& sys)
			{
				return sys;
			}


			/**
			Homogenize and patch the target system.
			*/
			static 
			void PrepareTarget(SystemType & target)
			{
				// target system came from the constructor
				target.Homogenize(); // work over projective coordinates
				target.AutoPatch(); // then patch if needed
			}

			static void FormStart(StartSystemType & start, SystemType const& target)
			{
				start = StartSystemType(target);	
			}

			static
			void FormHomotopy(SystemType & homotopy, SystemType const& target, StartSystemType const& start, std::string const& path_variable_name)
			{
				auto t = MakeVariable(path_variable_name); 

				homotopy = (1-t)*target + MakeRational(node::Rational::Rand())*t*start;
				homotopy.AddPathVariable(t);
			}

			/**
			\brief Sets up the homotopy for the system to be solved.
			
			1. Homogenizes the system, 
			2. patches it, 
			3. constructs the start system,
			4. stores the number of start points, 
			5. makes a path variable,
			6. forms the straight line homotopy between target and start, with the gamma trick

			boom, you're ready to go.
			*/
			void SystemSetup(std::string const& path_variable_name)
			{
				PrepareTarget(TargetSystem());

				// now we populate the start system
				FormStart(StartSystem(), TargetSystem());

				FormHomotopy(Homotopy(), TargetSystem(), StartSystem(), path_variable_name);
			}


			/**
			A getter for the system to be tracked to.
			*/
			SystemT& TargetSystem() 
			{
				return target_system_;
			}

			/**
			A getter for the homotopy being used.
			*/
			SystemT& Homotopy() 
			{
				return homotopy_;
			}

			/**
			A getter for the start system being used.
			*/
			StartSystemT & StartSystem() 
			{
				return start_system_;
			}
		};


		/**
		This system management policy allows the user to prevent the zero dim algorithm from making clones, and instead the burden of supplying the target system, start system, and homotopy are entirely up to the user.

		Using this policy implies the user manages these things entirely.

		\see CloneGiven
		*/
		template<typename SystemType, typename StartSystemType>
		struct RefToGiven : public SysMgmtPolicy<RefToGiven<SystemType, StartSystemType>, SystemType, StartSystemType, 
								std::reference_wrapper< const SystemType>, std::reference_wrapper< const StartSystemType>>
		{
			using SMP = SysMgmtPolicy<RefToGiven<SystemType, StartSystemType>, SystemType, StartSystemType, 
								std::reference_wrapper< const SystemType>, std::reference_wrapper< const StartSystemType>>;
			friend SMP;

			using StoredSystemT = typename SMP::StoredSystemT;
			using StoredStartSystemT = typename SMP::StoredStartSystemT;

			using SMP::TargetSystem;
			using SMP::StartSystem;
			using SMP::Homotopy;


private:
			StoredSystemT target_system_; ///< The target system which we track to.
			StoredStartSystemT start_system_; ///< The start system, which produces start points.
			StoredSystemT homotopy_; ///< homotopy, on which we wish the path vanishes.


public:
			/**
			Simply forward references to the given systems on the stored systems.
			*/
			RefToGiven(SystemType const& target, StartSystemType const& start, SystemType const& hom)
			 : 
			 	target_system_(std::ref(target)), 
			 	start_system_(std::ref(start)), 
			 	homotopy_(std::ref(hom))
			{}


			/**
			\brief Store a reference to the argument system.
			*/
			static
			StoredSystemT AtConstruct(SystemType const& sys)
			{
				return std::ref(sys);
			}

			/**
			\brief Store a reference to the argument system.
			*/
			template<typename T>
			static
			T AtSet(T const& sys)
			{
				return std::ref(sys);
			}

			void SystemSetup(std::string const& path_variable_name) const
			{ }

		};
	} // ns policy
} // ns bertini