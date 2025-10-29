#include "TerrainMesh.hpp"

TerrainMesh::TerrainMesh() 
    : vao_(0)
    , vbo_(0)
    , ebo_(0)
    , indexCount_(0) 
{

}

TerrainMesh::~TerrainMesh()
{
    cleanup();
}


void TerrainMesh::generateFromServer(TerrainServerOpenTopo& server,
                        double latStart, double lonStart,
                        double latEnd, double lonEnd,
                        int rows, int cols)
{
    std::cout << "Fetching terrain data...\n";
    auto points = server.getElevationGrid(latStart, lonStart, latEnd, lonEnd, rows, cols);
    
    // Find elevation range
    double minElev = 1e9, maxElev = -1e9;
    for (const auto& p : points)
    {
        minElev = std::min(minElev, p.alt);
        maxElev = std::max(maxElev, p.alt);
    }

    double elevRange = maxElev - minElev;
    std::cout << "Elevation range: " << minElev << " to " << maxElev << " m\n";

    // Generate vertices
    std::vector<Vertex> vertices;
    vertices.reserve(rows * cols);

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            int idx = i * cols + j;
            const auto& pt = points[idx];
            
            float x = static_cast<float>(j);
            float z = static_cast<float>(i);
            float y = (elevRange > 0) ? 
                        static_cast<float>((pt.alt - minElev) / elevRange) * 10.0f : 0.0f;

            // Height-based coloring
            float normHeight = (elevRange > 0) ? (pt.alt - minElev) / elevRange : 0.5f;
            glm::vec3 color = getTerrainColor(normHeight);

            vertices.push_back({
                glm::vec3(x, y, z),
                glm::vec3(0.0f, 1.0f, 0.0f), // Will calculate normals later
                color
            });
        }
    }

    // Generate indices for triangles
    std::vector<unsigned int> indices;
    indices.reserve((rows - 1) * (cols - 1) * 6);

    for (int i = 0; i < rows - 1; ++i)
    {
        for (int j = 0; j < cols - 1; ++j)
        {
            unsigned int topLeft = i * cols + j;
            unsigned int topRight = topLeft + 1;
            unsigned int bottomLeft = (i + 1) * cols + j;
            unsigned int bottomRight = bottomLeft + 1;

            // First triangle
            indices.push_back(topLeft);
            indices.push_back(bottomLeft);
            indices.push_back(topRight);

            // Second triangle
            indices.push_back(topRight);
            indices.push_back(bottomLeft);
            indices.push_back(bottomRight);
        }
    }

    // Calculate normals
    calculateNormals(vertices, indices, cols);

    // Create OpenGL buffers
    createBuffers(vertices, indices);

    // Center the mesh
    centerPoint_ = glm::vec3(cols / 2.0f, 5.0f, rows / 2.0f);
    
    std::cout << "Generated terrain mesh: " << vertices.size() << " vertices, "
                << indices.size() / 3 << " triangles\n";
}

void TerrainMesh::render(bool wireframe)
{
    if (vao_ == 0) return;

    glBindVertexArray(vao_);
    
    if (wireframe)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    }
    else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    }

    glDrawElements(GL_TRIANGLES, indexCount_, GL_UNSIGNED_INT, 0);
    
    glBindVertexArray(0);
}

void TerrainMesh::calculateNormals(std::vector<Vertex>& vertices,
                        const std::vector<unsigned int>& indices,
                        int cols)
{
    // Reset normals
    for (auto& v : vertices)
    {
        v.normal = glm::vec3(0.0f);
    }

    // Accumulate face normals
    for (size_t i = 0; i < indices.size(); i += 3)
    {
        unsigned int i0 = indices[i];
        unsigned int i1 = indices[i + 1];
        unsigned int i2 = indices[i + 2];

        glm::vec3 v0 = vertices[i0].position;
        glm::vec3 v1 = vertices[i1].position;
        glm::vec3 v2 = vertices[i2].position;

        glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

        vertices[i0].normal += normal;
        vertices[i1].normal += normal;
        vertices[i2].normal += normal;
    }

    // Normalize
    for (auto& v : vertices)
    {
        v.normal = glm::normalize(v.normal);
    }
}

glm::vec3 TerrainMesh::getTerrainColor(float height)
{
    // Realistic terrain coloring
    if (height < 0.2f)
        return glm::mix(glm::vec3(0.2f, 0.5f, 0.3f), glm::vec3(0.3f, 0.6f, 0.3f), height / 0.2f); // Dark to light green
    else if (height < 0.5f)
        return glm::mix(glm::vec3(0.3f, 0.6f, 0.3f), glm::vec3(0.6f, 0.5f, 0.3f), (height - 0.2f) / 0.3f); // Green to brown
    else if (height < 0.7f)
        return glm::mix(glm::vec3(0.6f, 0.5f, 0.3f), glm::vec3(0.5f, 0.5f, 0.5f), (height - 0.5f) / 0.2f); // Brown to gray
    else
        return glm::mix(glm::vec3(0.5f, 0.5f, 0.5f), glm::vec3(1.0f, 1.0f, 1.0f), (height - 0.7f) / 0.3f); // Gray to white (snow)
}

void TerrainMesh::createBuffers(const std::vector<Vertex>& vertices,
                    const std::vector<unsigned int>& indices)
{
    cleanup(); // Clean up old buffers if any

    indexCount_ = indices.size();

    // Generate and bind VAO
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // Generate and bind VBO
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex),
                vertices.data(), GL_STATIC_DRAW);

    // Generate and bind EBO
    glGenBuffers(1, &ebo_);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo_);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int),
                indices.data(), GL_STATIC_DRAW);

    // Position attribute
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, position));

    // Normal attribute
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, normal));

    // Color attribute
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex),
                        (void*)offsetof(Vertex, color));

    glBindVertexArray(0);
}

void TerrainMesh::cleanup()
{
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (ebo_) glDeleteBuffers(1, &ebo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    vbo_ = ebo_ = vao_ = 0;
}