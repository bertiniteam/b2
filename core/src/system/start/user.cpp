//This file is part of Bertini 2.
//
//mhom.cpp is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//mhom.cpp is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with mhom.cpp.  If not, see <http://www.gnu.org/licenses/>.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license,
// as well as COPYING.  Bertini2 is provided with permitted
// additional terms in the b2/licenses/ directory.

#include "bertini2/system/start/user.hpp"
#include "bertini2/detail/sha256.hpp"

#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>


BOOST_CLASS_EXPORT(bertini::start_system::User);


namespace {

    // a double as its exact IEEE-754 bit pattern, the convention of every b2 digest preimage
    void EmitExact(std::ostream& out, double v)
    {
        std::uint64_t bits;
        static_assert(sizeof bits == sizeof v, "a double is 64 bits");
        std::memcpy(&bits, &v, sizeof bits);
        out << "d64:" << std::hex << std::setw(16) << std::setfill('0') << bits << std::dec;
    }

    // an mpfr value at its own precision with every stored digit, as the system canonical
    // encoding writes its multiprecision constants
    template<typename RealT>
    void EmitExactMp(std::ostream& out, RealT const& v)
    {
        out << v.str(0, std::ios::scientific);
    }

} // namespace


namespace bertini {

    namespace start_system {

        // constructor for User start system, from any other *suitable* system.
        User::User(System const& s, SampCont<complex_dbl> const& solns) : user_system_(s), solns_in_dbl_(true)
        {
            std::get<SampCont<complex_dbl>>(solns_) = solns;
        }

        User::User(System const& s, SampCont<complex_mp> const& solns) : user_system_(s), solns_in_dbl_(false)
        {
            std::get<SampCont<complex_mp>>(solns_) = solns;
        }


        unsigned long long User::NumStartPoints() const
        {
            if (solns_in_dbl_)
                return std::get<SampCont<complex_dbl>>(solns_).size();
            else
                return std::get<SampCont<complex_mp>>(solns_).size();
        }


        std::string User::GivenStartIdentity() const
        {
            // the exact start data, in order: version, kind, point count, then each point's
            // length and coordinates.  Versioned so a later change of encoding is a new identity.
            std::ostringstream text;
            text << "b2start/1 ";
            if (solns_in_dbl_)
            {
                auto const& pts = std::get<SampCont<complex_dbl>>(solns_);
                text << "dbl " << pts.size() << '\n';
                for (auto const& p : pts)
                {
                    text << p.size();
                    for (Eigen::Index ii = 0; ii < p.size(); ++ii)
                    {
                        text << ' ';
                        EmitExact(text, p(ii).real());
                        text << ' ';
                        EmitExact(text, p(ii).imag());
                    }
                    text << '\n';
                }
            }
            else
            {
                auto const& pts = std::get<SampCont<complex_mp>>(solns_);
                text << "mp " << pts.size() << '\n';
                for (auto const& p : pts)
                {
                    text << p.size();
                    for (Eigen::Index ii = 0; ii < p.size(); ++ii)
                    {
                        text << ' ' << p(ii).precision() << ' ';
                        EmitExactMp(text, p(ii).real());
                        text << ' ';
                        EmitExactMp(text, p(ii).imag());
                    }
                    text << '\n';
                }
            }
            return detail::Sha256(text.str()).Hex();
        }



        Vec<complex_dbl> User::GenerateStartPoint(complex_dbl,unsigned long long index) const
        {
            if (solns_in_dbl_)
                return std::get<SampCont<complex_dbl>>(solns_)[index];
            else
            {
                const auto& r = std::get<SampCont<complex_mp>>(solns_)[index];
                Vec<complex_dbl> pt(r.size());
                for (unsigned ii=0; ii<r.size(); ++ii)
                    pt(ii) = complex_dbl(r(ii));

                return pt;
            }
        }


        Vec<complex_mp> User::GenerateStartPoint(complex_mp,unsigned long long index) const
        {
            if (solns_in_dbl_)
            {
                const auto& r = std::get<SampCont<complex_dbl>>(solns_)[index];
                Vec<complex_mp> pt(r.size());
                for (unsigned ii=0; ii<r.size(); ++ii)
                    pt(ii) = static_cast<complex_mp>(r(ii));

                return pt;
            }
            else
            {
                return std::get<SampCont<complex_mp>>(solns_)[index];
            }
        }

    } // namespace start_system
} //namespace bertini
