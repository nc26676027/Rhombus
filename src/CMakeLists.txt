set(RHOMBUS_SRC_LIST
        "matvec.cpp"
        "seal_api.cpp"
        "matrix.cpp"
        "matmul.cpp"
        "CtSerialization.cpp"
        )

add_library(Rhombus ${RHOMBUS_SRC_LIST})
target_link_libraries(Rhombus PUBLIC SEAL::seal)

# Bootstrap sub-module (RBOOT)
add_subdirectory(bootstrap)
