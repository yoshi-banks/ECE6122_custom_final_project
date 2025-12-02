/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file Camera.cpp
 * 
 * @brief Camera class implementation
 * @details This file implements the Camera class for handling 3D camera movements
 * and view matrix calculations in the Terrain Visualizer application.
 */

#include "Camera.hpp"

#include <GLFW/glfw3.h>

/**
 * @brief Constructor for Camera class
 */
Camera::Camera(glm::vec3 position,
        glm::vec3 target,
        glm::vec3 up)
    : position_(position)
    , target_(target)
    , worldUp_(up)
    , distance_(glm::length(position - target))
    , yaw_(-90.0f)
    , pitch_(0.0f)
    , mouseSensitivity_(0.3f)
    , zoomSensitivity_(5.0f)
    , movementSpeed_(50.0f)
    , panSpeed_(0.5f)
{
    updateCameraVectors();
}

/**
 * @brief Get the view matrix for the camera
 */
glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(position_, target_, up_);
}

/**
 * @brief Process mouse movement input
 */
void Camera::processMouseMovement(float xoffset, float yoffset, bool constrainPitch)
{
    xoffset *= mouseSensitivity_;
    yoffset *= mouseSensitivity_;

    yaw_ += xoffset;
    pitch_ += yoffset;

    if (constrainPitch)
    {
        if (pitch_ > 89.0f) pitch_ = 89.0f;
        if (pitch_ < -89.0f) pitch_ = -89.0f;
    }

    updateCameraVectors();
}

/**
 * @brief Process mouse scroll input
 */
void Camera::processMouseScroll(float yoffset)
{
    distance_ -= yoffset * zoomSensitivity_;
    if (distance_ < 1.0f) distance_ = 5.0f;         // Minimum zoom closer
    if (distance_ > 100.0f) distance_ = 500.0f;     // Maximum zoom further
    updateCameraVectors();
}

/**
 * @brief Process keyboard input for camera movement
 */
void Camera::processKeyboard(int key, float deltaTime)
{
    float velocity = movementSpeed_ * deltaTime;
    
    // Calculate movement directions relative to camera's current view
    // Use the actual front vector (not flattened) for true camera-relative movement
    glm::vec3 forward = front_;  // Direction camera is looking
    glm::vec3 right = right_;     // Camera's right vector
    glm::vec3 up = up_;           // Camera's up vector
    
    // WASD movement (moves the target in camera-relative directions)
    if (key == GLFW_KEY_W)
        target_ += forward * velocity;  // Move forward in look direction
    if (key == GLFW_KEY_S)
        target_ -= forward * velocity;  // Move backward from look direction
    if (key == GLFW_KEY_A)
        target_ -= right * velocity;    // Strafe left
    if (key == GLFW_KEY_D)
        target_ += right * velocity;    // Strafe right
    
    // QE for vertical movement (world up/down)
    if (key == GLFW_KEY_Q)
        target_ += up * velocity;
    if (key == GLFW_KEY_E)
        target_ -= up * velocity;
    
    // Arrow keys for rotation
    if (key == GLFW_KEY_LEFT)
        yaw_ -= 100.0f * deltaTime;
    if (key == GLFW_KEY_RIGHT)
        yaw_ += 100.0f * deltaTime;
    if (key == GLFW_KEY_UP)
        pitch_ += 100.0f * deltaTime;
    if (key == GLFW_KEY_DOWN)
        pitch_ -= 100.0f * deltaTime;
    
    // Constrain pitch
    if (pitch_ > 89.0f) pitch_ = 89.0f;
    if (pitch_ < -89.0f) pitch_ = -89.0f;
    
    // Zoom with Z/X
    if (key == GLFW_KEY_Z)
        distance_ -= velocity * 2.0f;
    if (key == GLFW_KEY_X)
        distance_ += velocity * 2.0f;
    
    distance_ = glm::clamp(distance_, 5.0f, 500.0f);
    
    updateCameraVectors();
}

/**
 * @brief Pan the camera
 */
void Camera::pan(float dx, float dy)
{
    glm::vec3 right = glm::normalize(glm::cross(front_, up_));
    glm::vec3 up = up_;
    
    target_ += right * dx * 0.1f;
    target_ += up * dy * 0.1f;
    updateCameraVectors();
}

/**
 * @brief Set the target of the camera
 */
void Camera::setTarget(glm::vec3 target)
{
    target_ = target;
    updateCameraVectors();
}

/**
 * @brief Set the distance from the camera to the target
 */
void Camera::setDistance(float distance)
{
    distance_ = distance;
    updateCameraVectors();
}

/**
 * @brief Set the movement speed of the camera
 */
void Camera::setMovementSpeed(float speed)
{
    movementSpeed_ = speed;
}

/**
 * @brief Set the pan speed of the camera
 */
void Camera::setPanSpeed(float speed)
{
    panSpeed_ = speed;
}

/**
 * @brief Set the rotation speed of the camera
 */
void Camera::setRotationSpeed(float speed)
{
    mouseSensitivity_ = speed;
}

/**
 * @brief Update the camera's vectors based on current yaw, pitch, and distance
 */
void Camera::updateCameraVectors()
{
    // Calculate position based on spherical coordinates around target
    glm::vec3 direction;
    direction.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    direction.y = sin(glm::radians(pitch_));
    direction.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    
    front_ = glm::normalize(direction);
    position_ = target_ - front_ * distance_;
    
    right_ = glm::normalize(glm::cross(front_, worldUp_));
    up_ = glm::normalize(glm::cross(right_, front_));
}