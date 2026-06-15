# Face module (encrypted 1:N face embedding recognition)
# Gated on RHOMBUS_BUILD_FACE option

# Phase 1: HyDia engine (slot-domain EE-MVM)
# set(FACE_SOURCES
#     hydia.cpp
#     diagonal_layout.cpp
# )

# Phase 2: Threshold + cleaning
# list(APPEND FACE_SOURCES
#     threshold.cpp
#     clean.cpp
# )

# Phase 3: Tag select + count
# list(APPEND FACE_SOURCES
#     tagselect.cpp
#     count.cpp
# )

# Placeholder: no sources yet, will be populated in Phase 1
set(FACE_SOURCES "")

if(FACE_SOURCES)
    add_library(RhombusFace STATIC ${FACE_SOURCES})
    target_include_directories(RhombusFace PUBLIC
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/src/face
    )
    target_link_libraries(RhombusFace PUBLIC Rhombus SEAL::seal)
endif()
