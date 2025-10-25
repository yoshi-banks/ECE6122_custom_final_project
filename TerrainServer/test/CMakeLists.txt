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

# Optional: Add tests to CTest
enable_testing()
add_test(NAME BasicTerrainTest COMMAND test_terrain_basic)