#include <Arduino.h>
#include "Vertex.h"

Vertex::Vertex() {
    pos = {0,0,0};
    norm = {0,0,0};
    color = {1,1,1};
}
Vertex::Vertex(Matrix<3> &_pos, Matrix<3> &_norm) {
    pos = _pos;
    norm = _norm;
}

Vertex::Vertex(Matrix<3> &_pos, Matrix<3> &_norm, Matrix<3> &_pos_s) {
    pos = _pos;
    norm = _norm;
    pos_s = _pos_s;
}

Vertex::Vertex(Matrix<3> &_pos, Matrix<3> &_norm, Matrix<3> &_pos_s, Matrix<3> &_color) {
    pos = _pos;
    norm = _norm;
    pos_s = _pos_s;
    color = _color;
}

Vertex::Vertex(float *_pos, float *_norm, float *_color) {
    pos = {_pos[0], _pos[1], _pos[2]};
    norm = {_norm[0], _norm[1], _norm[2]};
    color = {_color[0], _color[1], _color[2]};
}

Vertex::Vertex(const float *_pos, const float *_norm) {
    pos = {_pos[0], _pos[1], _pos[2]};
    norm = {_norm[0], _norm[1], _norm[2]};
}

bool Vertex::operator==(Vertex &other) {
    return pos(0) == other.pos(0) && pos(1) == other.pos(1) && pos(2) == other.pos(2);
}