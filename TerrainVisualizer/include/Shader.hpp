/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Shader.hpp
 * 
 * @brief Shader class definition
 * @details This file declares the Shader class for handling OpenGL shaders in the Terrain Visualizer application.
 */

#pragma once

#include <GL/glew.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

/**
 * @brief Shader class
 * 
 * @details This class represents an OpenGL shader program for rendering.
 */
class Shader
{
public:
    Shader();

    bool loadFromString(const std::string& vertexCode, const std::string& fragmentCode);
    void use() const;

    void setMat4(const std::string& name, const glm::mat4& mat) const;
    void setVec3(const std::string& name, const glm::vec3& value) const;

    void setFloat(const std::string& name, float value) const;
    void setInt(const std::string& name, int value) const;

    static Shader createDefaultTerrainShader();
    static Shader createCubeShader();

private:
    GLuint compileShader(const char* code, GLenum type);

    GLuint programID_;
};