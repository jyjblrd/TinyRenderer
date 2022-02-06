#ifndef Renderer_h
#define Renderer_h

#include <Arduino.h>
#include "Vertex.h"
#include "TFT_Driver.h"
#include "TTGO_T_Display.h"
#include "dragon.h"

#define MAX_VERTICES 300

#define RES_REDUX 1

#define Ka 0.1
#define Kd 0.5
#define Ks 0.7
#define SHININESS 20

#define RGB565
#define Z_BUFF_16

class Renderer {
    private:
        TFT_Driver tft;

        int num_vertices;
        int num_indices;

        Vertex shaded_vertices[dragon::vertex_count];

        #ifdef Z_BUFF_16
            uint16_t z_buffer[BLOCK_HEIGHT][BLOCK_WIDTH];
        #else
            uint8_t z_buffer[BLOCK_HEIGHT][BLOCK_WIDTH];
        #endif


        //uint16_t cur_frame_buffer[TFT_HEIGHT/RES_REDUX + 1][TFT_WIDTH/RES_REDUX + 1];
        #ifdef RGB565
            //uint16_t cur_frame_buffer[TFT_HEIGHT][TFT_WIDTH];
            uint16_t next_frame_buffer[BLOCK_HEIGHT][BLOCK_WIDTH];
        #else
            //uint8_t cur_frame_buffer[TFT_HEIGHT][TFT_WIDTH];
            uint8_t next_frame_buffer[BLOCK_HEIGHT][BLOCK_WIDTH];
        #endif

        //uint8_t should_erase[TFT_HEIGHT][TFT_WIDTH/8 + 1];
        //uint8_t pixel_changed[TFT_HEIGHT][TFT_WIDTH/8 + 1];

        Matrix<4> cam_eye; // Pos of camera
        Matrix<4> cam_look; // Center of screen
        Matrix<4> cam_up; // Up vector of camera

        Matrix<4> light_pos;
        Matrix<3> light_color;
        Matrix<3> trans_light_pos;

        float n = 0.09; // Near
        float f = 0.20; // Far
        float aspect = (float)TFT_WIDTH/TFT_HEIGHT;
        float fov = 2*PI * (30/360.0);

        float dist = 1; // Distance to screen center from camera

        Matrix<4,4> view_matrix;
        Matrix<4,4> proj_matrix;
        Matrix<4,4> model_matrix = {
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1
        };
         Matrix<4,4> view_model_matrix = view_matrix * model_matrix;
        Matrix<4,4> view_model_normal_matrix =~((view_model_matrix).Inverse());

    public:
        Renderer(TFT_Driver _tft);

        void set_vertices(const float *_vertices, const float *_normals, uint16_t size);
        void set_indices(const uint16_t *_indices, uint16_t size);

        Vertex vertex_shader(const float *vertex, const float *normal, Matrix<3> &color);
        void fragment_shader(uint32_t x, uint32_t y, uint32_t x_offset, uint32_t y_offset, Vertex *triag);

        void set_model_matrix(Matrix<4,4> &_model_matrix);

        void calc_view_matrix();

        void render();
        void draw_buffer();

        bool is_backface(Vertex *triag);

        float mag(Matrix<4> &m);
        Matrix<4> norm_homo(Matrix<4> &m);
        Matrix<3> from_homo(Matrix<4> const &m);
        float det(Matrix<2> u, Matrix<2> v);
        void triag_bbox(Vertex *triag, uint32_t &min_x, uint32_t &min_y, uint32_t &max_x, uint32_t &max_y);

        void baycentric_coords(Vertex *triag, Matrix<2> v, float &alpha, float &beta, float &gamma);
        float ilf(Matrix<2> a, Matrix<2> b, Matrix<2> v);

        Matrix<3> mult(float s, Matrix<3> const &m);
        Matrix<3> norm(Matrix<3> const &m);
        Matrix<3> cross(Matrix<3> const &a, Matrix<3> const &b);

        void task1(void * parameter);
};

#endif