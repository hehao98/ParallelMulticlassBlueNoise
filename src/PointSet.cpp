//
// Created by 何昊 on 2018/12/17.
//

#include "PointSet.h"

#include <random>
#include <iostream>

PointSet::PointSet()
{
    vbuf = nullptr;
    tbuf = nullptr;
    vao = vbo1 = vbo2 = 0;
    vbufSize = tbufSize = 0;
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
        Point p;
        p.pos.x = distX(generator);
        p.pos.y = distY(generator);
        p.type = 0;
        points.push_back(p);
    }
}

void PointSet::readPointSetFromFile(const char *path,
                                    double xmin, double xmax, double ymin, double ymax)
{
    std::ifstream infile;
    infile.open(path);
    if (infile.bad()) {
        std::cerr << "Bad Filepath " << path << std::endl;
        return;
    }

    this->xmin = xmin;
    this->xmax = xmax;
    this->ymin = ymin;
    this->ymax = ymax;
    points.clear();
    while (infile.good()) {
        int type;
        double x, y;
        infile >> type >> x >> y;
        Point p;
        p.pos = glm::dvec2(x, y);
        p.type = type;
        points.push_back(p);
    }
}

void PointSet::updateRenderData(float screenXmin, float screenXmax,
                                float screenYmin, float screenYmax)
{
    // Reconstruct vbuf to represent new data
    delete[] vbuf;
    delete[] tbuf;
    vbufSize = 2 * points.size();
    vbuf = new float[vbufSize];
    for (int i = 0; i < points.size(); ++i) {
        vbuf[2 * i]     = screenXmin + (screenXmax - screenXmin) *
                          (float)((points[i].pos.x - this->xmin)/(this->xmax - this->xmin));
        vbuf[2 * i + 1] = screenYmin + (screenYmax - screenYmin) *
                          (float)((points[i].pos.y - this->ymin)/(this->ymax - this->ymin));
    }
    tbufSize = points.size();
    tbuf = new int[tbufSize];
    for (int i = 0; i < tbufSize; ++i) {
        tbuf[i] = points[i].type;
    }

    // Updata data in GPU
    if (vao == 0 && vbo1 == 0 && vbo2 == 0) {
        glGenBuffers(1, &vbo1);
        glGenBuffers(1, &vbo2);
        glGenVertexArrays(1, &vao);

        glBindVertexArray(vao);
            glBindBuffer(GL_ARRAY_BUFFER, vbo1);
                glBufferData(GL_ARRAY_BUFFER, vbufSize * sizeof(float), vbuf, GL_DYNAMIC_DRAW);
                glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), nullptr);
                glEnableVertexAttribArray(0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);

            glBindBuffer(GL_ARRAY_BUFFER, vbo2);
                glBufferData(GL_ARRAY_BUFFER, tbufSize * sizeof(int), tbuf, GL_DYNAMIC_DRAW);
                glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(int), nullptr);
                glEnableVertexAttribArray(1);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    } else {
        glBindBuffer(GL_ARRAY_BUFFER, vbo1);
        glBufferSubData(GL_ARRAY_BUFFER, 0, vbufSize * sizeof(float), vbuf);
        glBindBuffer(GL_ARRAY_BUFFER, vbo2);
        glBufferSubData(GL_ARRAY_BUFFER, 0, tbufSize * sizeof(int), tbuf);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void PointSet::render(const Shader &shader)
{
    shader.use();
    glBindVertexArray(vao);
    glDrawArrays(GL_POINTS, 0, vbufSize / 2);
}


