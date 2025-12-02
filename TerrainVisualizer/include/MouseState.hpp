/**
 * Author: Joshua Miller
 * Class: ECE6122 (Q)
 * Last Date Modified: 2025-12-01
 * 
 * @file MouseState.hpp
 * @brief MouseState struct definition
 */

#pragma once

struct MouseState
{
    bool firstMouse = true;
    bool leftButtonPressed = false;
    bool rightButtonPressed = false;
    bool middleButtonPressed = false;
    double lastX = 400.0;
    double lastY = 300.0;
};