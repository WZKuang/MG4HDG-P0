if(NETGEN_USE_PYTHON)
    add_ngsolve_python_module(ngsprol_py 
      python_ngsprol.cpp
      ../python/python_ngsprol.cpp
      ../utils/NCprolongation.cpp
      ../utils/facetfespace0.cpp
    )
    set_target_properties(ngsprol_py PROPERTIES INSTALL_RPATH "${NETGEN_RPATH_TOKEN}/../${NETGEN_PYTHON_RPATH}")
    install(TARGETS ngsprol_py DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/prol COMPONENT prol)

    target_include_directories(ngsprol_py PUBLIC ${NGSOLVE_INCLUDE_DIRS})
    install (FILES
        __init__.py 
        DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/prol
        COMPONENT prol
    )

# build stub files for pybind11 packages
if(BUILD_STUB_FILES)
find_program(PYBIND11_STUBS NAMES pybind11-stubgen)
if(PYBIND11_STUBS)
  message("-- Found pybind11-stubgen: ${PYBIND11_STUBS}")
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/prol-stubs/)
  install(CODE "execute_process(COMMAND ${PYBIND11_STUBS} --no-setup-py prol)")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../stubs/prol-stubs/ DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/prol/)
else(PYBIND11_STUBS)
  message(WARNING "pybind11-stubgen not found, if you want to create stub files
for better autocompletion support install it with pip.")
endif(PYBIND11_STUBS)
endif(BUILD_STUB_FILES)

endif(NETGEN_USE_PYTHON)
