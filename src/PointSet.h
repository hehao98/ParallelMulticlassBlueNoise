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
    explicit PointSet();

    void generateWhiteNoisePointSet(int size,
                                    double xmin, double xmax, double ymin, double ymax);

    void updateRenderData(float screenXmin, float screenXmax,
                          float screenYmin, float screenYmax);

    void render(const Shader &shader);
private:
    std::vector<glm::dvec2> points;
    double xmin, xmax, ymin, ymax;

    // Members for OpenGL draw calls
    float *vbuf;
    int bufSize;
    GLuint vao, vbo;
};


#endif //PROJECT_POINTSET_H
