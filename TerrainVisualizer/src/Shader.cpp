/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Shader.cpp
 * @brief Shader class implementation
 * @details This file implements the Shader class for handling OpenGL shaders
 * in the Terrain Visualizer application.
 */

#include "Shader.hpp"

/**
 * @brief Constructor for Shader class
 */
Shader::Shader()
    : programID_(0) 
{

}

/**
 * @brief Loads a shader program from strings
 * @param vertexCode Vertex shader code
 * @param fragmentCode Fragment shader code
 */
bool Shader::loadFromString(const std::string& vertexCode, const std::string& fragmentCode)
{
    // Compile vertex shader
    GLuint vertexShader = compileShader(vertexCode.c_str(), GL_VERTEX_SHADER);
    if (vertexShader == 0) return false;

    // Compile fragment shader
    GLuint fragmentShader = compileShader(fragmentCode.c_str(), GL_FRAGMENT_SHADER);
    if (fragmentShader == 0)
    {
        glDeleteShader(vertexShader);
        return false;
    }

    // Link shaders
    programID_ = glCreateProgram();
    glAttachShader(programID_, vertexShader);
    glAttachShader(programID_, fragmentShader);
    glLinkProgram(programID_);

    // Check for linking errors
    GLint success;
    glGetProgramiv(programID_, GL_LINK_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetProgramInfoLog(programID_, 512, nullptr, infoLog);
        std::cerr << "Shader linking failed:\n" << infoLog << "\n";
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        return false;
    }

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    return true;
}

/**
 * @brief Uses the shader program
 */
void Shader::use() const
{
    glUseProgram(programID_);
}

/**
 * @brief Sets a mat4 uniform
 * @param name Name of the uniform
 * @param mat Matrix value
 */
void Shader::setMat4(const std::string& name, const glm::mat4& mat) const
{
    glUniformMatrix4fv(glGetUniformLocation(programID_, name.c_str()),
                        1, GL_FALSE, glm::value_ptr(mat));
}

/**
 * @brief Sets a vec3 uniform
 * @param name Name of the uniform
 * @param value Vector value
 */
void Shader::setVec3(const std::string& name, const glm::vec3& value) const
{
    glUniform3fv(glGetUniformLocation(programID_, name.c_str()),
                1, glm::value_ptr(value));
}

/**
 * @brief Sets a float uniform
 * @param name Name of the uniform
 * @param value Float value
 */
void Shader::setFloat(const std::string& name, float value) const
{
    glUniform1f(glGetUniformLocation(programID_, name.c_str()), value);
}

/**
 * @brief Sets an int uniform
 * @param name Name of the uniform
 * @param value Int value
 */
void Shader::setInt(const std::string& name, int value) const
{
    glUniform1i(glGetUniformLocation(programID_, name.c_str()), value);
}

/**
 * @brief Compiles a shader
 * @param code Shader source code
 * @param type Shader type (GL_VERTEX_SHADER or GL_FRAGMENT_SHADER)
 * @return Compiled shader ID, or 0 on failure
 */
GLuint Shader::compileShader(const char* code, GLenum type)
{
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &code, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr << "Shader compilation failed (" 
                    << (type == GL_VERTEX_SHADER ? "vertex" : "fragment")
                    << "):\n" << infoLog << "\n";
        return 0;
    }

    return shader;
}

/**
 * @brief Creates a default terrain shader
 */
Shader Shader::createDefaultTerrainShader()
{
    const std::string vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;
        layout (location = 2) in vec3 aColor;

        out vec3 FragPos;
        out vec3 Normal;
        out vec3 Color;

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;

        void main()
        {
            FragPos = vec3(model * vec4(aPos, 1.0));
            Normal = mat3(transpose(inverse(model))) * aNormal;
            Color = aColor;
            gl_Position = projection * view * model * vec4(aPos, 1.0);
        }
    )";

    const std::string fragmentShaderSource = R"(
        #version 330 core
        out vec4 FragColor;

        in vec3 FragPos;
        in vec3 Normal;
        in vec3 Color;

        uniform vec3 lightPos;
        uniform vec3 viewPos;
        uniform vec3 lightColor;

        void main()
        {
            // Ambient
            float ambientStrength = 0.3;
            vec3 ambient = ambientStrength * lightColor;

            // Diffuse
            vec3 norm = normalize(Normal);
            vec3 lightDir = normalize(lightPos - FragPos);
            float diff = max(dot(norm, lightDir), 0.0);
            vec3 diffuse = diff * lightColor;

            // Specular
            float specularStrength = 0.2;
            vec3 viewDir = normalize(viewPos - FragPos);
            vec3 reflectDir = reflect(-lightDir, norm);
            float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
            vec3 specular = specularStrength * spec * lightColor;

            // Green enhancement for flat/top surfaces
            // How much the normal points upward (Y-axis)
            float upwardness = max(dot(norm, vec3(0.0, 1.0, 0.0)), 0.0);
            upwardness = pow(upwardness, 0.5); // Soften the effect

            // Add green tint to flatter surfaces
            vec3 greenTint = vec3(0.2, 0.4, 0.2) * upwardness * 0.6;
            vec3 enhancedColor = Color + greenTint;

            vec3 result = (ambient + diffuse + specular) * enhancedColor;
            FragColor = vec4(result, 1.0);
        }
    )";

    Shader shader;
    shader.loadFromString(vertexShaderSource, fragmentShaderSource);
    return shader;
}

/**
 * @brief Creates a simple cube shader
 */
Shader Shader::createCubeShader()
{
    const std::string vertexShaderSource = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aColor;

        out vec3 FragColor;

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;

        void main()
        {
            FragColor = aColor;
            gl_Position = projection * view * model * vec4(aPos, 1.0);
        }
    )";

    const std::string fragmentShaderSource = R"(
        #version 330 core
        out vec4 FragColorOut;
        in vec3 FragColor;

        void main()
        {
            FragColorOut = vec4(FragColor, 1.0);
        }
    )";

    Shader shader;
    shader.loadFromString(vertexShaderSource, fragmentShaderSource);
    return shader;
}