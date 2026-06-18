set(BOOTSTRAP_SRC_LIST
    "Bootstrapper.cpp"
    "ModularReducer.cpp"
    "Polynomial.cpp"
    "PolyUpdate.cpp"
    "Remez.cpp"
    "RemezApp.cpp"
    "MinicompFunc.cpp"
    "MinicompRemez.cpp"
    "Choosemax.cpp"
    "Point.cpp"
    "func.cpp"
    "comp/program.cpp"
    "comp/SEALfunc.cpp"
)

add_library(RhombusBootstrap ${BOOTSTRAP_SRC_LIST})
target_include_directories(RhombusBootstrap PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Find NTL library (high-precision arithmetic, required by RBOOT bootstrapping)
find_package(NTL QUIET)
if(NOT NTL_FOUND)
    # Fallback: search homebrew paths on macOS
    find_path(NTL_INCLUDE_DIR NTL/RR.h
        PATHS /opt/homebrew/include /usr/local/include)
    find_library(NTL_LIBRARY ntl
        PATHS /opt/homebrew/lib /usr/local/lib)
    if(NTL_INCLUDE_DIR AND NTL_LIBRARY)
        set(NTL_FOUND TRUE)
        add_library(NTL::ntl UNKNOWN IMPORTED)
        set_target_properties(NTL::ntl PROPERTIES
            IMPORTED_LOCATION "${NTL_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${NTL_INCLUDE_DIR}")
    endif()
endif()

if(NOT NTL_FOUND)
    message(FATAL_ERROR "NTL library not found. Install with: brew install ntl")
endif()

target_link_libraries(RhombusBootstrap PUBLIC SEAL::seal NTL::ntl)

# Data generation tool: produces arcsin_d127.txt needed by bootstrapping
add_executable(gen_arcsin gen_arcsin.cpp)
target_link_libraries(gen_arcsin PRIVATE RhombusBootstrap SEAL::seal NTL::ntl)

