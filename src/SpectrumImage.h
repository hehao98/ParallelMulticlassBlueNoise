//
// Created by 何昊 on 2019/1/16.
//

#ifndef PARALLELMULTICLASSBLUENOISE_SPECTRUMIMAGE_H
#define PARALLELMULTICLASSBLUENOISE_SPECTRUMIMAGE_H

#include <glad/glad.h>
#include <vector>

#include "Shader.h"
#include "PointSet.h"

class SpectrumImage {
public:
    void updateSpectrum(std::vector<PointSet::Point> points);
    void render(const Shader &shader);
private:
    const static int size = 1024;
    float data[size][size];
    GLuint texture;
    GLuint VAO, VBO, EBO;
};

#endif //PARALLELMULTICLASSBLUENOISE_SPECTRUMIMAGE_H
