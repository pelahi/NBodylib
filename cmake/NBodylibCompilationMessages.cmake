#
# Tell the world what what we are doing
#
macro(nbody_report feature)

	# Output feature name and underscore it in the next line
	message("\n${feature}")
	string(REGEX REPLACE "." "-" _underscores ${feature})
	message("${_underscores}\n")

	set(_args "${ARGN}")
	list(LENGTH _args _nargs)
	math(EXPR _nargs "${_nargs} - 1")
	foreach(_idx RANGE 0 ${_nargs} 2)

		# Items in the list come with a message first, then the variable name
		list(GET _args ${_idx} _msg)
		math(EXPR _idx2 "${_idx} + 1")
		list(GET _args ${_idx2} _varname)

		# We try to keep things up to 80 cols
		string(LENGTH ${_msg} _len)
		math(EXPR _nspaces "75 - ${_len}")
		string(RANDOM LENGTH ${_nspaces} _spaces)
		string(REGEX REPLACE "." " " _spaces "${_spaces}")
		string(CONCAT _msg "${_msg}" ${_spaces})
		message(" ${_msg} ${NBODY_HAS_${_varname}}")
	endforeach()
endmacro()

macro(nbody_compilation_summary)
	message("\nNBodyLib successfully configured with the following settings:")
	nbody_report(Dependencies "OpenMP" OPENMP
	    "OpenMP target directive for GPU accelaration supported and active" OPENMPGPU)
	nbody_report(Dependencies "OpenACC" OPENACC)
	nbody_report(Types "All calculations/properties stored as float" SINGLE_PRECISION "All integeres are long int" LONG_INT)
	nbody_report("Particle data" "Do not store mass, all particles are the same mass" NO_MASS
	             "Use single precision to store positions, velocities, other props" SINGLE_PARTICLE_PRECISION
	             "Use unsigned particle PIDs" UNSIGNED_PARTICLE_PIDS "Use unsigned particle IDs" UNSIGNED_PARTICLE_IDS
	             "Activate gas" USE_GAS "Activate stars" USE_STARS "Activate black holes/sinks" USE_BH
	             "Activate extra dm properties" USE_EXTRA_DM_PROPERTIES
	             "Extra input info stored" USE_EXTRA_INPUT_INFO
	             "Extra FOF info stored" USE_EXTRA_FOF_INFO
	             "Large memory KDTree" USE_LARGE_KDTREE
	             "Particle compiled for SWIFT" USE_SWIFT_INTERFACE
			 )
	message("")
	message("Compilation")
	message("=========================")
	message("Compiler ID          : ${CMAKE_CXX_COMPILER_ID}" )
	message("Compiler             : ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_COMPILER_ID}" )
	message("Build Type           : ${CMAKE_BUILD_TYPE}")
	message("Build Flags          : ${CMAKE_CXX_${CMAKE_BUILD_TYPE}_FLAGS}")
	message("----------------------")
	message(" Include directories : ${NBODY_INCLUDE_DIRS}")
	message(" Macros defined      : ${NBODY_DEFINES}")
	message(" Libs                : ${NBODY_LIBS}")
	message(" C++ flags           : ${NBODY_CXX_FLAGS}")
	message(" Link flags          : ${NBODY_LINK_FLAGS}")
	message("")
endmacro()
