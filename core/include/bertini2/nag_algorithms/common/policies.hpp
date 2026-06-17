//This file is part of Bertini 2.
//
//bertini2/nag_algorithms/common/policies.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/nag_algorithms/common/policies.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/nag_algorithms/common/policies.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, 
// as well as COPYING.  Bertini2 is provided with permitted 
// additional terms in the b2/licenses/ directory.

// individual authors of this file include:
// silviana amethyst, university of wisconsin-eau claire

/**
\file bertini2/nag_algorithms/common/policies.hpp 

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

			// This policy owns (deep-copies) its systems, so in a distributed solve rank 0's
			// systems can be broadcast and installed authoritatively on every rank.  RefToGiven,
			// which only holds references to user-managed systems, sets this false.
			static constexpr bool OwnsSystems = true;

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
				// MakeHomotopy builds H = (1-t)*target + gamma*t*start with a random gamma,
				// choosing a blend block when the start system carries structured blocks (e.g.
				// the MHom products-of-linears start) and node arithmetic otherwise.  The same
				// construction is exposed to Python as system.make_homotopy so a user-authored
				// start system can be turned into a trackable homotopy.
				homotopy = MakeHomotopy(target, start, path_variable_name);
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

			// The user owns these systems (we only hold references); a distributed solve must not
			// overwrite them, so it does not broadcast/install rank 0's systems here.  See CloneGiven.
			static constexpr bool OwnsSystems = false;


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

			void SystemSetup(std::string const& /*path_variable_name*/) const
			{ }

		};




		/**
		\brief A base class for system management for algorithms which operate on a single system.

		In contrast to SysMgmtPolicy, some algorithms -- notably the NumericalIrreducibleDecomposition algorithm -- do not form a homotopy from a start system, and so only need to manage a single target system.  This base provides just the target-system get/set.

		\see SysMgmtPolicy
		*/
		template< typename D, typename SystemType, typename StoredSystemType>
		struct SingleSysMgmtPolicy
		{
			using SystemT = SystemType;
			using StoredSystemT = StoredSystemType;

private:
			const D& AsDerived() const
			{
				return static_cast<const D&>(*this);
			}

			D& AsDerived()
			{
				return static_cast<D&>(*this);
			}

public:
			/**
			A getter for the system to be operated on.
			*/
			const SystemT& TargetSystem() const
			{
				return AsDerived().target_system_;
			}

			/**
			A setter for the system to be operated on.
			*/
			void TargetSystem(StoredSystemT const& sys)
			{
				AsDerived().target_system_ = sys;
			}
		};




		/**
		\brief A single-system management policy which clones the supplied target.

		The single-system analog of CloneGiven: makes a clone of the supplied system at construct time, so the algorithm may homogenize/patch it without disturbing the user's original.  No start system or homotopy is formed.

		Reusable by any algorithm which operates on a single system with no homotopy-from-start-system (e.g. NumericalIrreducibleDecomposition).

		\see CloneGiven, SingleSysMgmtPolicy
		*/
		template<typename SystemType>
		struct CloneTarget : public SingleSysMgmtPolicy<CloneTarget<SystemType>, SystemType, SystemType>
		{
			using SMP = SingleSysMgmtPolicy<CloneTarget<SystemType>, SystemType, SystemType>;
			friend SMP;

			using StoredSystemT = typename SMP::StoredSystemT;

			using SMP::TargetSystem;

			using SystemT = SystemType;

private:
			StoredSystemT target_system_;

public:
			/**
			Forward the system on to the stored (cloned) target.
			*/
			CloneTarget(SystemType const& target) : target_system_(AtConstruct(target))
			{}


			static
			StoredSystemT AtConstruct(SystemType const& sys)
			{
				return Clone(sys);
			}


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
				target.Homogenize(); // work over projective coordinates
				target.AutoPatch(); // then patch if needed
			}


			/**
			\brief Sets up the (single) system to be operated on.

			Homogenizes and patches the target.  No start system or homotopy is formed.
			*/
			void SystemSetup()
			{
				PrepareTarget(TargetSystem());
			}


			/**
			A non-const getter for the system to be operated on.
			*/
			SystemT& TargetSystem()
			{
				return target_system_;
			}
		};
	} // ns policy
} // ns bertini