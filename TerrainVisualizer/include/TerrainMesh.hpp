#pragma once

#include <vector>
#include <GL/glew.h>
#include <glm/glm.hpp>
#include <iostream>
#include "TerrainServerOpenTopo.hpp"

struct Vertex
{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec3 color;
};

class TerrainMesh
{
public:
    TerrainMesh();
    
    ~TerrainMesh();

    void generateFromServer(TerrainServerOpenTopo& server,
                           double latStart, double lonStart,
                           double latEnd, double lonEnd,
                           int rows, int cols);

    void render(bool wireframe = false);

    glm::vec3 getCenterPoint() const { return centerPoint_; }

private:
    void calculateNormals(std::vector<Vertex>& vertices,
                         const std::vector<unsigned int>& indices,
                         int cols);

    glm::vec3 getTerrainColor(float height);

    void createBuffers(const std::vector<Vertex>& vertices,
                      const std::vector<unsigned int>& indices);

    void cleanup();

    GLuint vao_;
    GLuint vbo_;
    GLuint ebo_;
    unsigned int indexCount_;
    glm::vec3 centerPoint_;
};