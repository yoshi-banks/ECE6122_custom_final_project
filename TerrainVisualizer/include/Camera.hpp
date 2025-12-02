/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Camera.hpp
 * @brief Camera class definition
 * 
 * @details This file declares the Camera class for handling 3D camera movements
 * and view matrix calculations in the Terrain Visualizer application.
 */

#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

/**
 * @brief Camera class
 * 
 * @details This class represents a 3D camera for handling camera movements and view matrix calculations.
 */
class Camera
{
public:
    Camera(glm::vec3 position = glm::vec3(0.0f, 10.0f, 20.0f),
           glm::vec3 target = glm::vec3(0.0f, 0.0f, 0.0f),
           glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f));

    glm::mat4 getViewMatrix() const;

    void processMouseMovement(float xoffset, float yoffset, bool constrainPitch = true);
    void processKeyboard(int key, float deltaTime);
    void processMouseScroll(float yoffset);
    void setTarget(glm::vec3 target);
    void setDistance(float distance);
    void pan(float dx, float dy);

    glm::vec3 getPosition() const { return position_; }
    glm::vec3 getTarget() const { return target_; }
    glm::vec3 getFront() const { return front_; }

    void setMovementSpeed(float speed);
    void setPanSpeed(float speed);
    void setRotationSpeed(float speed);

private:
    void updateCameraVectors();

    glm::vec3 position_;
    glm::vec3 target_;
    glm::vec3 front_;
    glm::vec3 up_;
    glm::vec3 right_;
    glm::vec3 worldUp_;

    float distance_;
    float yaw_;
    float pitch_;

    float mouseSensitivity_;
    float zoomSensitivity_;
    float movementSpeed_;
    float panSpeed_;
};