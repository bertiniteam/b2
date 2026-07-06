//This file is part of Bertini 2.
//
//bertini2/parallel/mpi_utils.hpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//bertini2/parallel/mpi_utils.hpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with bertini2/parallel/mpi_utils.hpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team

/**
\file bertini2/parallel/mpi_utils.hpp

\brief Serialization-aware MPI send/recv helpers using plain C MPI + Boost.Serialization.

Boost.MPI's send/recv automatically invoked Boost.Serialization for complex types.
These helpers replicate that behavior using the C MPI API directly, removing the
Boost::mpi component dependency while keeping all existing serialize() methods.

Pattern: serialize to Boost binary_oarchive → single MPI_Send(MPI_BYTE).
         MPI_Probe + MPI_Get_count → MPI_Recv(MPI_BYTE) → deserialize.

The probe-before-receive approach avoids the two-message (size + data) race condition
that would arise with MPI_ANY_SOURCE.
*/

#pragma once

#ifdef BERTINI2_HAVE_MPI

#include "bertini2/parallel/mpi_include.hpp"

#include <sstream>
#include <string>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace bertini {
namespace parallel {


/**
\brief Serialize \p obj via Boost binary archive and send as a single MPI_BYTE message.
*/
template<typename T>
void mpi_send_serialized(MPI_Comm comm, int dest, int tag, T const& obj)
{
	std::ostringstream oss;
	{
		boost::archive::binary_oarchive oa(oss);
		oa << obj;
	}
	std::string const buf = oss.str();
	MPI_Send(buf.data(), static_cast<int>(buf.size()), MPI_BYTE, dest, tag, comm);
}


/**
\brief Receive a serialized message from any source, deserialize into \p obj.

Uses MPI_Probe to determine the source and byte count before receiving, so
this is safe to call with MPI_ANY_SOURCE even when the message is a single send.

\return The rank of the sending process.
*/
template<typename T>
int mpi_recv_serialized_any(MPI_Comm comm, int tag, T& obj)
{
	MPI_Status probe_status;
	MPI_Probe(MPI_ANY_SOURCE, tag, comm, &probe_status);

	int source = probe_status.MPI_SOURCE;
	int count  = 0;
	MPI_Get_count(&probe_status, MPI_BYTE, &count);

	std::string buf(static_cast<std::size_t>(count), '\0');
	MPI_Status recv_status;
	MPI_Recv(buf.data(), count, MPI_BYTE, source, tag, comm, &recv_status);

	std::istringstream iss(buf);
	boost::archive::binary_iarchive ia(iss);
	ia >> obj;

	return source;
}


/**
\brief Receive a serialized message from a specific source, deserialize into \p obj.
*/
template<typename T>
void mpi_recv_serialized(MPI_Comm comm, int source, int tag, T& obj)
{
	MPI_Status probe_status;
	MPI_Probe(source, tag, comm, &probe_status);

	int count = 0;
	MPI_Get_count(&probe_status, MPI_BYTE, &count);

	std::string buf(static_cast<std::size_t>(count), '\0');
	MPI_Status recv_status;
	MPI_Recv(buf.data(), count, MPI_BYTE, source, tag, comm, &recv_status);

	std::istringstream iss(buf);
	boost::archive::binary_iarchive ia(iss);
	ia >> obj;
}


/**
\brief Broadcast a std::string from \p root to all ranks in \p comm.

Sends length then content using plain MPI — no Boost involved.
On non-root ranks, \p s is overwritten with the broadcast value.
*/
inline void mpi_broadcast_string(MPI_Comm comm, std::string& s, int root)
{
	int len = static_cast<int>(s.size());
	MPI_Bcast(&len, 1, MPI_INT, root, comm);
	s.resize(static_cast<std::size_t>(len));
	MPI_Bcast(s.data(), len, MPI_CHAR, root, comm);
}


/**
\brief Broadcast a Boost-serializable object from \p root to all ranks in \p comm.

\p root serializes \p obj and broadcasts the bytes; every other rank deserializes into \p obj
(overwriting it).  Used to make \p root the single authoritative source of an object (e.g. the
homotopy / start system) rather than having each rank re-derive its own copy.
*/
template<typename T>
void mpi_broadcast_serialized(MPI_Comm comm, T& obj, int root)
{
	int rank = 0;
	MPI_Comm_rank(comm, &rank);

	std::string s;
	if (rank == root)
	{
		std::ostringstream oss;
		boost::archive::binary_oarchive oa(oss);
		oa << obj;
		s = oss.str();
	}

	mpi_broadcast_string(comm, s, root);

	if (rank != root)
	{
		std::istringstream iss(s);
		boost::archive::binary_iarchive ia(iss);
		ia >> obj;
	}
}


} // namespace parallel
} // namespace bertini

#endif // BERTINI2_HAVE_MPI
