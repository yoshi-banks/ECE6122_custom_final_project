/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Cube.cpp
 * 
 * @brief Cube class implementation
 * 
 * @details This file implements the Cube class for rendering 3D cubes in OpenGL.
 * Cubes can be created using latitude/longitude (converted to NED internally)
 * or NED coordinates directly.
 */

#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

#include "Cube.hpp"

/**
 * @brief Constructor for Cube class using lat/lon
 * 
 * @param mesh TerrainMesh object
 * @param lat Latitude in degrees
 * @param lon Longitude in degrees
 * @param color Color of the cube
 * @param size Size of the cube
 */
Cube::Cube(const TerrainMesh& mesh, double lat, double lon,
           const glm::vec3& color, float size)
    : size_(size), color_(color), vao_(0), vbo_(0), ebo_(0)
{
    // Get terrain elevation at this lat/lon
    double terrain_alt = mesh.getElevationAtLatLon(lat, lon);
    
    // Convert lat/lon/alt to NED
    NEDPoint ned = mesh.latLonToNED(lat, lon, terrain_alt);
    
    // Position in rendering coordinates (North->X, -Down->Y, East->Z)
    position_ = glm::vec3(
        static_cast<float>(ned.north),
        static_cast<float>(-ned.down),  // Negative because down is positive downward
        static_cast<float>(ned.east)
    );
    
    createBuffers();
}

/**
 * @brief Constructor for Cube class using NED coordinates
 * 
 * @param position NEDPoint object
 * @param color Color of the cube
 * @param size Size of the cube
 */
Cube::Cube(const NEDPoint& position, const glm::vec3& color, float size)
    : size_(size), color_(color), vao_(0), vbo_(0), ebo_(0)
{
    // Convert NED to rendering coordinates
    position_ = glm::vec3(
        static_cast<float>(position.north),
        static_cast<float>(-position.down),
        static_cast<float>(position.east)
    );
    
    createBuffers();
}

/**
 * @brief Destructor for Cube class
 * 
 * @details Deletes the OpenGL buffers associated with the cube.
 */
Cube::~Cube()
{
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (ebo_) glDeleteBuffers(1, &ebo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
}

/**
 * @brief Move constructor
 * @param other The other Cube to move from
 */
Cube::Cube(Cube&& other) noexcept
    : position_(other.position_)
    , color_(other.color_)
    , size_(other.size_)
    , vao_(other.vao_)
    , vbo_(other.vbo_)
    , ebo_(other.ebo_)
{
    // Prevent other from deleting the OpenGL buffers
    other.vao_ = 0;
    other.vbo_ = 0;
    other.ebo_ = 0;
}

/**
 * @brief Move assignment operator
 * @param other The other Cube to move from
 */
Cube& Cube::operator=(Cube&& other) noexcept
{
    if (this != &other) {
        // Clean up our existing OpenGL resources
        if (vbo_) glDeleteBuffers(1, &vbo_);
        if (ebo_) glDeleteBuffers(1, &ebo_);
        if (vao_) glDeleteVertexArrays(1, &vao_);
        
        // Transfer ownership from other
        position_ = other.position_;
        color_ = other.color_;
        size_ = other.size_;
        vao_ = other.vao_;
        vbo_ = other.vbo_;
        ebo_ = other.ebo_;
        
        // Prevent other from deleting the resources we now own
        other.vao_ = 0;
        other.vbo_ = 0;
        other.ebo_ = 0;
    }
    return *this;
}

/**
 * @brief Create OpenGL buffers for the cube
 */
void Cube::createBuffers()
{
    // Cube vertices (8 corners)
    float s = size_ / 2.0f;
    std::vector<Vertex> vertices = {
        // Front face (z+)
        {{-s, -s,  s}, color_},  // 0
        {{ s, -s,  s}, color_},  // 1
        {{ s,  s,  s}, color_},  // 2
        {{-s,  s,  s}, color_},  // 3
        // Back face (z-)
        {{-s, -s, -s}, color_},  // 4
        {{ s, -s, -s}, color_},  // 5
        {{ s,  s, -s}, color_},  // 6
        {{-s,  s, -s}, color_},  // 7
    };

    // Indices for 12 triangles (6 faces)
    std::vector<unsigned int> indices = {
        // Front
        0, 1, 2,  2, 3, 0,
        // Right
        1, 5, 6,  6, 2, 1,
        // Back
        5, 4, 7,  7, 6, 5,
        // Left
        4, 0, 3,  3, 7, 4,
        // Top
        3, 2, 6,  6, 7, 3,
        // Bottom
        4, 5, 1,  1, 0, 4
    };

    // Create VAO
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // Create VBO
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), 
                vertices.data(), GL_STATIC_DRAW);

    // Create EBO
    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
                indices.data(), GL_STATIC_DRAW);

    // Position attribute (layout location 0)
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), 
                         (void*)offsetof(Vertex, position));

    // Color attribute (layout location 1)  
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                         (void*)offsetof(Vertex, color));

    glBindVertexArray(0);
}

/**
 * @brief Render the cube
 * 
 * @param model_matrix The model matrix to apply 
 */
void Cube::render(const glm::mat4& model_matrix)
{
    if (vao_ == 0) return;

    // Apply local transform (translate to cube position)
    glm::mat4 local_model = glm::translate(model_matrix, position_);
    
    // Set model matrix uniform (assumes shader is already active)
    GLint shader_program;
    glGetIntegerv(GL_CURRENT_PROGRAM, &shader_program);
    GLint modelLoc = glGetUniformLocation(shader_program, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, &local_model[0][0]);

    glBindVertexArray(vao_);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}