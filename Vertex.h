#ifndef Vertex_h
#define Vertex_h

#include <Arduino.h>
#include <BasicLinearAlgebra.h>
using namespace BLA;

class Vertex {
    public:
        Vertex();
        Vertex(Matrix<3> &_pos, Matrix<3> &_norm);
        Vertex(Matrix<3> &_pos, Matrix<3> &_norm, Matrix<3> &_pos_s);
        Vertex(Matrix<3> &_pos, Matrix<3> &_norm, Matrix<3> &_pos_s, Matrix<3> &_color);
        Vertex(float *_pos, float *_norm, float *_color);
        Vertex(const float *_pos, const float *_norm);

        Matrix<3> pos; // View space
        Matrix<3> pos_s; // Position on screen coordinate space
        Matrix<3> norm; // View space
        Matrix<3> color = {1,1,1};

        bool operator == (Vertex &other);
};

#endif