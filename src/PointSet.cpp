//
// Created by 何昊 on 2018/12/17.
//

#include "PointSet.h"

#include <random>

PointSet::PointSet()
{
    vbuf = nullptr;
    vao = vbo = 0;
    bufSize = 0;
}

void PointSet::generateWhiteNoisePointSet(int size,
                                          double xmin, double xmax, double ymin, double ymax)
{
    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distX(xmin, xmax);
    std::uniform_real_distribution<double> distY(ymin, ymax);

    points.clear();
    for (int i = 0; i < size; ++i) {
        glm::dvec2 v;
        v.x = distX(generator);
        v.y = distY(generator);
        points.push_back(v);
    }
}

void PointSet::updateRenderData(float screenXmin, float screenXmax,
                                float screenYmin, float screenYmax)
{
    // Reconstruct vbuf to represent new data
    if (vbuf != nullptr) delete[] vbuf;
    bufSize = 2 * points.size();
    vbuf = new float[bufSize];
    for (int i = 0; i < points.size(); ++i) {
        vbuf[2 * i]     = screenXmin + (screenXmax - screenXmin) *
                          (float)((points[i].x - this->xmin)/(this->xmax - this->xmin));
        vbuf[2 * i + 1] = screenYmin + (screenYmax - screenYmin) *
                          (float)((points[i].y - this->ymin)/(this->ymax - this->ymin));
    }

    // Updata data in GPU
    if (vao == 0 && vbo == 0) {
        glGenBuffers(1, &vao);
        glGenVertexArrays(1, &vbo);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferData(GL_ARRAY_BUFFER, bufSize * sizeof(float), vbuf, GL_DYNAMIC_DRAW);
            glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
            glEnableVertexAttribArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    } else {
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
            glBufferSubData(GL_ARRAY_BUFFER, 0, bufSize * sizeof(float), vbuf);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void PointSet::render(const Shader &shader)
{
    shader.use();
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, bufSize / 2);
}


