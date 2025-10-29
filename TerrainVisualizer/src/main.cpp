#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>

#include "TerrainServerOpenTopo.hpp"
#include "Camera.hpp"
#include "TerrainMesh.hpp"
#include "Shader.hpp"

// Global state for mouse input
struct MouseState
{
    bool firstMouse = true;
    bool leftButtonPressed = false;
    bool rightButtonPressed = false;
    double lastX = 400.0;
    double lastY = 300.0;
} mouseState;

Camera camera;

// Callbacks
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (mouseState.firstMouse)
    {
        mouseState.lastX = xpos;
        mouseState.lastY = ypos;
        mouseState.firstMouse = false;
    }

    double xoffset = xpos - mouseState.lastX;
    double yoffset = mouseState.lastY - ypos; // Reversed since y-coordinates go from bottom to top

    mouseState.lastX = xpos;
    mouseState.lastY = ypos;

    if (mouseState.leftButtonPressed)
    {
        camera.processMouseMovement(xoffset, yoffset);
    }
    else if (mouseState.rightButtonPressed)
    {
        camera.pan(-xoffset, yoffset);
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT)
    {
        mouseState.leftButtonPressed = (action == GLFW_PRESS);
    }
    else if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
        mouseState.rightButtonPressed = (action == GLFW_PRESS);
    }
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.processMouseScroll(yoffset);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }
}

void printHelp()
{
    std::cout << "\n=== Terrain Visualizer Controls ===\n";
    std::cout << "Left Mouse + Drag:  Rotate camera\n";
    std::cout << "Right Mouse + Drag: Pan camera\n";
    std::cout << "Mouse Wheel:        Zoom in/out\n";
    std::cout << "ESC:                Exit\n";
    std::cout << "====================================\n\n";
}

int main(int argc, char** argv)
{
    printHelp();

    // Initialize GLFW
    if (!glfwInit())
    {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    // Set OpenGL version (3.3 Core)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create window
    const int WINDOW_WIDTH = 1200;
    const int WINDOW_HEIGHT = 800;
    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT,
                                          "Terrain Visualizer", nullptr, nullptr);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetKeyCallback(window, key_callback);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK)
    {
        std::cerr << "Failed to initialize GLEW\n";
        return -1;
    }

    // OpenGL configuration
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.53f, 0.81f, 0.92f, 1.0f); // Sky blue background

    // Initialize terrain server
    std::cout << "Initializing terrain server...\n";
    TerrainServerOpenTopo server("./terrain_cache", "https://api.opentopodata.org/v1", "srtm90m");
    server.setInterpolationTolerance(0.1);

    // Create terrain mesh
    TerrainMesh terrain;
    
    // Parse command line arguments or use defaults
    double latStart = 36.056361;
    double lonStart = -112.175607;
    double latEnd = 36.137089;
    double lonEnd = -112.105038;
    int rows = 50;
    int cols = 50;

    if (argc >= 5)
    {
        latStart = std::stod(argv[1]);
        lonStart = std::stod(argv[2]);
        latEnd = std::stod(argv[3]);
        lonEnd = std::stod(argv[4]);
    }
    if (argc >= 7)
    {
        rows = std::stoi(argv[5]);
        cols = std::stoi(argv[6]);
    }

    std::cout << "Loading terrain: [" << latStart << ", " << lonStart << "] to ["
              << latEnd << ", " << lonEnd << "] with " << rows << "x" << cols << " grid\n";

    terrain.generateFromServer(server, latStart, lonStart, latEnd, lonEnd, rows, cols);

    // Setup camera to look at terrain center
    camera.setTarget(terrain.getCenterPoint());

    // Create shader
    Shader shader = Shader::createDefaultTerrainShader();

    // Render loop
    std::cout << "\nStarting render loop...\n";
    
    double lastTime = glfwGetTime();
    int frameCount = 0;

    while (!glfwWindowShouldClose(window))
    {
        // Calculate FPS
        double currentTime = glfwGetTime();
        frameCount++;
        if (currentTime - lastTime >= 1.0)
        {
            std::cout << "FPS: " << frameCount << "\r" << std::flush;
            frameCount = 0;
            lastTime = currentTime;
        }

        // Clear buffers
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Use shader
        shader.use();

        // Set up transforms
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, -terrain.getCenterPoint()); // Center the terrain
        
        glm::mat4 view = camera.getViewMatrix();
        
        glm::mat4 projection = glm::perspective(
            glm::radians(45.0f),
            static_cast<float>(WINDOW_WIDTH) / static_cast<float>(WINDOW_HEIGHT),
            0.1f, 1000.0f
        );

        shader.setMat4("model", model);
        shader.setMat4("view", view);
        shader.setMat4("projection", projection);

        // Set lighting
        glm::vec3 lightPos = camera.getPosition() + glm::vec3(10.0f, 20.0f, 10.0f);
        shader.setVec3("lightPos", lightPos);
        shader.setVec3("viewPos", camera.getPosition());
        shader.setVec3("lightColor", glm::vec3(1.0f, 1.0f, 0.95f));

        // Render terrain
        terrain.render(false);

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Cleanup
    glfwTerminate();
    std::cout << "\nExiting...\n";
    
    return 0;
}