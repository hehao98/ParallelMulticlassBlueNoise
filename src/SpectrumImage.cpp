//
// Created by 何昊 on 2019/1/16.
//

#include "SpectrumImage.h"
#include <cmath>

void SpectrumImage::updateSpectrum(std::vector<PointSet::Point> points)
{
    if (texture == 0) {
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        // Set texture wrapping/filtering options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }

    float value = 0.001f;
    const float PI = 3.1415926f;
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < size; ++y) {
            data[x][y] = 0;
            for (int i = 0; i < points.size(); ++i) {
                data[x][y] += value * cos(2*PI*(points[i].pos.x*x + points[i].pos.y*y)/size);
            }
        }
    }

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, size, size, 0, GL_RED, GL_FLOAT, data);
    glGenerateMipmap(GL_TEXTURE_2D);
}

void SpectrumImage::readSpectrumFromFile(const char *filepath)
{
    std::ifstream infile;
    infile.open(filepath);
    if (infile.bad()) {
        std::cerr << "Bad Filepath " << filepath << std::endl;
        return;
    }

    std::string tmp;
    infile >> tmp;
    int dimension, height, width;
    infile >> dimension >>  height >>width;

    std::vector<float> image;
    for (int i = 0; i < height * width; ++i) {
        float tmp;
        infile >> tmp >> tmp >> tmp;
        image.push_back(tmp/255.0);
    }

    if (texture == 0) {
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        // Set texture wrapping/filtering options
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    }

    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, height, width, 0, GL_RED, GL_FLOAT, image.data());
    glGenerateMipmap(GL_TEXTURE_2D);
}

void SpectrumImage::render(const Shader &shader)
{
    if (VAO == 0) {
        float vertices[] = {
                // Positions        // colors         // texture coordinates
                0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f,
                0.0f, -1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,
                1.0f, -1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 0.0f,
                1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f,

        };
        unsigned int indices[] = {
                0, 1, 2, // first triangle
                0, 2, 3, // second triangle
        };

        // Vertex Buffer Objects(VBO) are used to pass vertices to GPU for the vertex shader
        // 1 is assigned as the unique ID to this VBO
        glGenBuffers(1, &VBO);

        // Vertex Array Objects(VAO) are used to store vertex attribute pointers
        glGenVertexArrays(1, &VAO);

        // Element Buffer Objects(EBO) describe additional information like indices of triangles
        glGenBuffers(1, &EBO);

        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);
        // Tell OpenGL how tp interpret vertex data
        // Pass data to layout(location=0), each data 3 values, type float, no normalization,
        // with the stride as 3*sizeof(float)
        // location 0: vertex location
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) 0);
        glEnableVertexAttribArray(0);
        // location 1: vertex color
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) (3 * sizeof(float)));
        glEnableVertexAttribArray(1);
        // location 2: texture coordinate
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void *) (6 * sizeof(float)));
        glEnableVertexAttribArray(2);
    }

    shader.use();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, texture);

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
}