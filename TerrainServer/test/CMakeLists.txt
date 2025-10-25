# Test executables for TerrainServer

# Test 1: Basic terrain data retrieval
add_executable(test_terrain_basic
    test_terrain_basic.cpp
)

target_link_libraries(test_terrain_basic PRIVATE
    terrain_server_lib
)

# Test 2: TCP Server test
add_executable(test_terrain_server
    test_terrain_server.cpp
)

target_link_libraries(test_terrain_server PRIVATE
    terrain_server_lib
)

# Test 3: Grid performance test
add_executable(test_terrain_grid
    test_terrain_grid.cpp
)

target_link_libraries(test_terrain_grid PRIVATE
    terrain_server_lib
)

# Test 4: Smart cache and interpolation test
add_executable(test_smart_cache
    test_smart_cache.cpp
)

target_link_libraries(test_smart_cache PRIVATE
    terrain_server_lib
)

# Test 5: High res Grand Canyon
add_executable(test_grand_canyon_high_res
    test_grand_canyon_high_res.cpp
)

target_link_libraries(test_grand_canyon_high_res PRIVATE
    terrain_server_lib
)

# Optional: Add tests to CTest
enable_testing()
add_test(NAME BasicTerrainTest COMMAND test_terrain_basic)
add_test(NAME SmartCacheTest COMMAND test_smart_cache)