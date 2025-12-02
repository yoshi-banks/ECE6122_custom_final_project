/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Cube.hpp
 * @brief Cube class definition
 * @details This file declares the Cube class for rendering cubes in the Terrain Visualizer application.
 * It uses OpenGL to draw the cube on the screen.
 */

#pragma once

#include <glm/glm.hpp>
#include <GL/glew.h>
#include <vector>

#include "TerrainMesh.hpp"

/**
 * @brief Cube class
 * @details This class represents a 3D cube for rendering in OpenGL.
 */
class Cube {
public:
    struct Vertex {
        glm::vec3 position;
        glm::vec3 color;
    };
    
    // Constructor using lat/lon (converts to NED internally)
    Cube(const TerrainMesh& mesh, double lat, double lon,
         const glm::vec3& color, float size = 5.0f);
    
    // Constructor using NED coordinates directly
    Cube(const NEDPoint& position, const glm::vec3& color, float size = 5.0f);
    
    ~Cube();

    // Move semantics (OpenGL handles shouldn't be copied)
    Cube(Cube&& other) noexcept;
    Cube& operator=(Cube&& other) noexcept;
    
    // Disable copying (OpenGL handles shouldn't be copied)
    Cube(const Cube&) = delete;
    Cube& operator=(const Cube&) = delete;
    
    void render(const glm::mat4& model_matrix);
    glm::vec3 getPosition() const { return position_; }

    void setPosition(const NEDPoint& pos) {
        position_ = glm::vec3(pos.north, -pos.down, pos.east);
    }
    
private:
    void createBuffers();
    
    glm::vec3 position_;  // Position in NED coordinates (north, -down, east)
    glm::vec3 color_;
    float size_;
    GLuint vao_;
    GLuint vbo_;
    GLuint ebo_;
};