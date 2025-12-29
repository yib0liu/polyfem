# Polyfem Solvers (https://github.com/polyfem/paraviewo)
# License: MIT

if(TARGET paraviewo::paraviewo)
    return()
endif()

message(STATUS "Third-party: creating target 'paraviewo::paraviewo'")

include(CPM)
CPMAddPackage("gh:yib0liu/paraviewo#84c48b2b6de3c28565cd3af83e3f14d9b3521e65")
