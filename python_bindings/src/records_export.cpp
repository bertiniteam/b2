//This file is part of Bertini 2.
//
//python_bindings/src/records_export.cpp is free software: you can redistribute it
//and/or modify it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//It is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
//See the GNU General Public License for more details.
//
// Copyright(C) Bertini2 Development Team
//
// See <http://www.gnu.org/licenses/> for a copy of the license, as well as COPYING.
// Bertini2 is provided with permitted additional terms in the b2/licenses/ directory.

/** \file records_export.cpp
\brief The `records` submodule: bindings for the structured output directory.

The directory format is plain text (that is the point: docs/records/b2rec-1.md),
so Python READS it with json/pandas, no bindings needed.  These bindings exist for
WRITING through the one C++ implementation -- session-journal appends (one writer per
file), content-addressed definitions, and the derived-view renderers -- so both faces
of bertini2 produce byte-identical structure.
*/

#include <boost/json.hpp>

#include "records_export.hpp"
#include <bertini2/records/output_directory.hpp>

namespace bertini { namespace python {

namespace {

	std::shared_ptr<records::OutputDirectory> MakeDirectory(std::string const& path)
	{
		// the process-shared instance: one session history file per directory per
		// process, however many OutputDirectory objects Python constructs
		return records::OutputDirectory::Shared(path);
	}

	void AppendJson(records::OutputDirectory& self, std::string const& record_json)
	{
		self.Append(boost::json::parse(record_json).as_object());
	}

	std::string PutDefinitionWrapper(records::OutputDirectory& self,
	                                 std::string const& content, std::string const& kind,
	                                 std::string const& external_id, std::string const& label)
	{
		return self.PutDefinition(content, kind,
		                          external_id.empty() ? std::nullopt
		                                              : std::optional<std::string>(external_id),
		                          label);
	}

	void AnnotateJson(records::OutputDirectory& self, std::string const& run_id,
	                  long long index, std::string const& key, std::string const& value_json)
	{
		self.Annotate(run_id, static_cast<std::int64_t>(index), key,
		              boost::json::parse(value_json));
	}

	std::string RootString(records::OutputDirectory const& self)
	{
		return self.Root().string();
	}

} // unnamed namespace

void ExportRecords()
{
	using namespace boost::python;

	scope current_scope;
	std::string new_submodule_name(extract<const char*>(current_scope.attr("__name__")));
	new_submodule_name.append(".records");
	object new_submodule(borrowed(PyImport_AddModule(new_submodule_name.c_str())));
	current_scope.attr("records") = new_submodule;

	scope new_submodule_scope = new_submodule;
	new_submodule_scope.attr("__doc__") =
		"The structured output directory (record schema b2rec/1): durable, plain-text "
		"records of computations.  READ it with json/pandas -- no bertini required; these "
		"bindings are for WRITING through the single C++ implementation.";

	class_<records::OutputDirectory, std::shared_ptr<records::OutputDirectory>, boost::noncopyable>
		("OutputDirectory", no_init)
		.def("__init__", make_constructor(&MakeDirectory),
			"Open (creating if needed) the structured output directory at the given path, "
			"writing its self-documenting README.txt on first creation.")
		.def("root", &RootString, (arg("self")), "The directory's root path.")
		.def("append", &AppendJson, (arg("self"), arg("record_json")),
			"Append one record (a JSON object as a string) to this session's history file.  "
			"One writer per file; flushed per record.")
		.def("put_definition", &PutDefinitionWrapper,
			(arg("self"), arg("content"), arg("kind"), arg("external_id") = "", arg("label") = ""),
			"Store a definition under the given kind folder ('systems'/'configs'/'givens'); "
			"returns the id (SHA-256 of the bytes unless external_id names it).  An optional "
			"label weaves a role into the filename (e.g. 'start_points') -- presentation only.")
		.def("annotate", &AnnotateJson,
			(arg("self"), arg("run"), arg("index"), arg("key"), arg("value_json")),
			"Attach metadata to a recorded point: appends an annotation record for "
			"({run, index}) with the given key and JSON-encoded value.  Newest wins "
			"per (point, key); annotations render into results.json beside the point.")
		.def("refresh_results", &records::OutputDirectory::RefreshResults, (arg("self")),
			"(Re)write the pretty-printed, self-complete results.json from the declared "
			"result records.")
		.def("refresh_index", &records::OutputDirectory::RefreshIndex, (arg("self")),
			"(Re)write INDEX.txt: one line per run.")
		.def("describe", &records::OutputDirectory::Describe, (arg("self")),
			"One human line: how much is here.")
		;
}

}} // namespaces
