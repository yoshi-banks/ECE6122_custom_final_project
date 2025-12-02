/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file PathLine.hpp
 * 
 * @brief PathLine class definition
 * @details This file declares the PathLine class for rendering a flight path line
 * in the Terrain Visualizer application using OpenGL.
 */

#pragma once

#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "NEDPoint.hpp"
#include "TerrainMesh.hpp"

/**
 * @brief PathLine class
 * 
 * @details This class represents a flight path line for rendering in OpenGL.
 */
class PathLine {
public:
    PathLine(const std::vector<NEDPoint>& path, const TerrainMesh& terrain) {
        // Convert NED path to vertex positions
        std::vector<glm::vec3> positions;
        for (const auto& ned : path) {
            float x = static_cast<float>(ned.north);
            float y = static_cast<float>(-ned.down);
            float z = static_cast<float>(ned.east);
            positions.push_back(glm::vec3(x, y, z));
        }
        
        // Create line VAO/VBO
        glGenVertexArrays(1, &vao_);
        glGenBuffers(1, &vbo_);
        
        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glBufferData(GL_ARRAY_BUFFER, positions.size() * sizeof(glm::vec3),
                    positions.data(), GL_STATIC_DRAW);
        
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
        
        glBindVertexArray(0);
        
        vertex_count_ = positions.size();
    }
    
    ~PathLine() {
        if (vbo_) glDeleteBuffers(1, &vbo_);
        if (vao_) glDeleteVertexArrays(1, &vao_);
    }
    
    void render() {
        if (vao_ == 0) return;
        
        glBindVertexArray(vao_);
        glLineWidth(3.0f);
        glDrawArrays(GL_LINE_STRIP, 0, vertex_count_);
        glBindVertexArray(0);
    }
    
private:
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    int vertex_count_ = 0;
};