//
// A class that describes a sampled point set
//
// Created by 何昊 on 2018/12/17.
//

#ifndef PROJECT_POINTSET_H
#define PROJECT_POINTSET_H

#include <vector>

// OpenGL Math Library
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glad/glad.h>

#include "Shader.h"

class PointSet {
public:
    struct Point {
        glm::dvec2 pos;
        int type; // start from zero
    };
    std::vector<Point> points;

    explicit PointSet();

    void generateWhiteNoisePointSet(int size,
                                    double xmin, double xmax, double ymin, double ymax);

    void readPointSetFromFile(const char *path,
                              double xmin, double xmax, double ymin, double ymax);

    void updateRenderData(float screenXmin, float screenXmax,
                          float screenYmin, float screenYmax, int renderClass = -1);

    void render(const Shader &shader);

private:
    double xmin, xmax, ymin, ymax;

    // Members for OpenGL draw calls
    float *vbuf; // vertex buffer
    int   *tbuf; // type buffer
    size_t vbufSize, tbufSize;
    GLuint vao, vbo1, vbo2;
};


#endif //PROJECT_POINTSET_H
