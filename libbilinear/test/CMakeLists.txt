set(bilinear_test_sources
  ElementSizesTest.cpp
  SerializeTest.cpp
  BasicViabilityTest.cpp
  TestFastMultExpDeRooij.cpp
)

foreach(appSrc ${bilinear_test_sources})
    get_filename_component(appName ${appSrc} NAME_WE)
    set(appDir ../bin/test)

    add_executable(${appName} ${appSrc})
    target_link_libraries(${appName} PRIVATE bilinear_mainapp bilinear)

    add_test(NAME ${appName} COMMAND ${appName})
    set_target_properties(${appName} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${appDir})
endforeach()
