set(SRC_FILES
        spring.cpp
        util.cpp
      )

add_library(tinyed SHARED ${SRC_FILES} ${INC_FILES} )

target_include_directories(tinyed PUBLIC
        ${PROJECT_SOURCE_DIR}/includes
        )

target_link_libraries(
    tinyed
    Eigen3::Eigen
)
