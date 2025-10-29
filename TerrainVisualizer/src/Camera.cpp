#include "Camera.hpp"

Camera::Camera(glm::vec3 position,
        glm::vec3 target,
        glm::vec3 up)
    : position_(position), target_(target), worldUp_(up),
        distance_(glm::length(position - target)),
        yaw_(-90.0f), pitch_(0.0f),
        mouseSensitivity_(0.1f), zoomSensitivity_(1.0f)
{
    updateCameraVectors();
}

glm::mat4 Camera::getViewMatrix() const
{
    return glm::lookAt(position_, target_, up_);
}

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

void Camera::processMouseScroll(float yoffset)
{
    distance_ -= yoffset * zoomSensitivity_;
    if (distance_ < 1.0f) distance_ = 1.0f;
    if (distance_ > 100.0f) distance_ = 100.0f;
    updateCameraVectors();
}

void Camera::setTarget(glm::vec3 target)
{
    target_ = target;
    updateCameraVectors();
}

void Camera::pan(float dx, float dy)
{
    glm::vec3 right = glm::normalize(glm::cross(front_, up_));
    glm::vec3 up = up_;
    
    target_ += right * dx * 0.1f;
    target_ += up * dy * 0.1f;
    updateCameraVectors();
}

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