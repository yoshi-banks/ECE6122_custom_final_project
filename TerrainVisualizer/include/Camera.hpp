#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

class Camera
{
public:
    Camera(glm::vec3 position = glm::vec3(0.0f, 10.0f, 20.0f),
           glm::vec3 target = glm::vec3(0.0f, 0.0f, 0.0f),
           glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f));

    glm::mat4 getViewMatrix() const;

    void processMouseMovement(float xoffset, float yoffset, bool constrainPitch = true);

    void processMouseScroll(float yoffset);

    void setTarget(glm::vec3 target);

    void pan(float dx, float dy);

    glm::vec3 getPosition() const { return position_; }
    glm::vec3 getTarget() const { return target_; }
    glm::vec3 getFront() const { return front_; }

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
};