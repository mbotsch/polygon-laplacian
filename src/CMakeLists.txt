set(HEADERS MeanCurvature.h  PolyDiffGeo.h GeodesicsInHeat.h  Parameterization.h Smoothing.h  SphericalHarmonics.h)

set(SOURCES MeanCurvature.cpp PolyDiffGeo.cpp GeodesicsInHeat.cpp  Parameterization.cpp Smoothing.cpp  SphericalHarmonics.cpp)

add_library(poly STATIC ${HEADERS} ${SOURCES})
target_link_libraries(poly pmp)

if(EMSCRIPTEN)
    add_executable(polylaplace main.cpp Viewer.cpp Viewer.h)
    target_link_libraries(polylaplace poly pmp_vis)
    set_target_properties(
        polylaplace
        PROPERTIES
          LINK_FLAGS
          "--shell-file ${PROJECT_SOURCE_DIR}/external/pmp-library/src/apps/data/shell.html --preload-file ${PROJECT_SOURCE_DIR}/data/bunny-hexa.off@input.off"
      )
else()
    add_executable(polylaplace main.cpp Viewer.cpp Viewer.h)
    target_link_libraries(polylaplace poly pmp_vis)

    add_executable(run_tests run_tests.cpp)
    target_link_libraries(run_tests poly pmp pmp_vis)
endif()