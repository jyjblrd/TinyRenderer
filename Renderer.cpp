#include <Arduino.h>
#include "Renderer.h"
#include "Vertex.h"
#include "TTGO_T_Display.h"

Renderer::Renderer(TFT_Driver _tft) {
    tft = _tft;

    for (size_t x = 0; x < BLOCK_WIDTH; x++) {
        for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
            z_buffer[y][x] = INFINITY;
            //cur_frame_buffer[y/RES_REDUX][x/RES_REDUX] = 0x00;
        }
    }

    cam_eye = {4, 0, 1, 1};
    cam_look = {0, 0, 0, 1};
    cam_up = {-2, 0, 1, 1};

    light_pos = {4,1,1,1};
    light_color = {1,1,1};

    float e = 1/tanf(fov/2);
    proj_matrix = {
        e/aspect,0,0,0,
        0,e,0,0,
        0,0,(f+n)/(n-f),(2*f*n)/(n-f),
        0,0,-1,0
    };
    // proj_matrix = {
    //     1/aspect,0,0,0,
    //     0,1,0,0,
    //     0,0,1,0,
    //     0,0,0,1,
    // };

    num_vertices = dragon::vertex_count;
    num_indices = dragon::vertex_count;
}

// void Renderer::set_vertices(const float *_vertices, const float *_normals, uint16_t size) {
//     if (size > MAX_VERTICES) {
//         Serial << "ERROR: too many verticies";    
//         return;
//     }

//     for (size_t i = 0; i < dragon::vertex_count; i++) {
//         Vertex v = Vertex(dragon::vertices[i], dragon::normals[i/3]);
//         vertices[i] = v;
//     }

//     num_vertices = size;
// }

// void Renderer::set_indices(const uint16_t *_indices, uint16_t size) {
//     if (size/3 > MAX_VERTICES) {
//         Serial << "ERROR: too many dragon::indices";    
//         return;
//     }

//     for (size_t i = 0; i < size; i++) {
//         dragon::indices[i] = _indices[i];
//     }   

//     num_indices = size;
// }

void Renderer::calc_view_matrix() {

    Matrix<4,4> t = {
        1,0,0,-cam_eye(0),
        0,1,0,-cam_eye(1),
        0,0,1,-cam_eye(2),
        0,0,0,1
    };

    Matrix<4> el = norm_homo(cam_look)-norm_homo(cam_eye);
    el(3) = 1; // Subtracting homogeneous cooridnates
    // Serial << "el" << el << "\n";
    float scale = dist/mag(el);
    Matrix<4,4> s = {
        scale,0,0,0,
        0,scale,0,0,
        0,0,scale,0,
        0,0,0,1
    };

    Matrix<4> lpp = s*t*cam_look; // l''
    float lppz = lpp(2)/lpp(3); // Remember that it is homogeneous coordinates
    float lppx = lpp(0)/lpp(3);
    float theta = -acosf( -lppz / sqrtf(lppx*lppx + lppz*lppz) );
    if (-lppz == sqrtf(lppx*lppx + lppz*lppz)) theta = -acosf(1); // For 0/0 case
    // Serial << "lpp " << lpp << "\n";
    // Serial << "lppz " << lppz << "\n";
    // Serial << "lppx " << lppx << "\n";
    // Serial << "theta " << theta << "\n";

    Matrix<4,4> r1 = {
        cosf(theta),0,sinf(theta),0,
        0,1,0,0,
        -sinf(theta),0,cosf(theta),0,
        0,0,0,1
    };

    Matrix<4> lppp = r1*lpp; // l'''
    float lpppz = lppp(2)/lppp(3);
    float lpppy = lppp(1)/lppp(3);
    float phi = -acosf( -lpppz / sqrtf(lpppy*lpppy + lpppz*lpppz) );
    if (-lpppz == sqrtf(lpppy*lpppy + lpppz*lpppz)) phi = -acosf(1); // For 0/0 case
    Matrix<4,4> r2 = {
        1,0,0,0,
        0,cosf(phi),-sinf(phi),0,
        0,sinf(phi),cosf(phi),0,
        0,0,0,1
    };

    Matrix<4> upppp = r2*r1*cam_up; // u'''' 
    float uppppy = upppp(1)/upppp(3);
    float uppppx = upppp(0)/upppp(3);
    float psi = -acosf( uppppy / sqrtf(uppppx*uppppx + uppppy*uppppy) );
    if (uppppy == sqrtf(uppppx*uppppx + uppppy*uppppy)) psi = -acosf(1); // For 0/0 case
    Matrix<4,4> r3 = {
        cosf(psi),-sinf(psi),0,0,
        sinf(psi),cosf(psi),0,0,
        0,0,1,0,
        0,0,0,1
    };

    view_matrix = r3*r2*r1*s*t;
    // Serial << r3 << "\n";
    // Serial << r2 << "\n";
    // Serial << r1 << "\n";
    // Serial << s << "\n";
    // Serial << t << "\n";
    // Serial << view_matrix << "\n";
}

void Renderer::set_model_matrix(Matrix<4,4> &_model_matrix) {
    model_matrix = _model_matrix;
}

Vertex Renderer::vertex_shader(const float *vertex, const float *normal, Matrix<3> &color) {

    // Position tranform
    Matrix<4> homo_pos = {vertex[0], vertex[1], vertex[2], 1}; // Convert to homogeneous coordinates

    homo_pos = view_model_matrix * homo_pos; // Project into view space
    Matrix<4> homo_pos_s = proj_matrix * homo_pos; // Project into projection space
    Matrix<3> pos = from_homo(homo_pos); // Homogeneous -> normal coords
    Matrix<3> pos_s = from_homo(homo_pos_s); // Homogeneous -> normal coords

    // Normal transform
    Matrix<4> homo_n = {normal[0], normal[1], normal[2], 1}; // Convert to homogeneous coordinates
    homo_n = view_model_normal_matrix * homo_n; // Project onto view coordinates
    Matrix<3> n = from_homo(homo_n); // Homogeneous -> normal coords

    // Phong reflection model to calculate ambiant and diffuse component
    Matrix<3> lm = norm(trans_light_pos - pos); // Light direction
    Matrix<3> rm = norm(-mult(2*(~lm*n)(0), n) + lm); // Reflected light
    Matrix<3> v = norm(-pos); // viewer

    Matrix<3> ambient = mult(Ka, color);
    Matrix<3> diffuse = mult(Kd*max((~lm*n)(0), (float)0), color);
    Matrix<3> specular = mult(Ks*pow(max((~rm*v)(0), (float)0), SHININESS), light_color);
    Matrix<3> vertex_color = ambient + diffuse + specular;

    Vertex res = Vertex(pos, n, pos_s, vertex_color);

    return res;
}

void Renderer::fragment_shader(uint32_t x, uint32_t y, uint32_t x_offset, uint32_t y_offset, Vertex *triag) {
    int start_time = micros();
    
    Matrix<2> v = {2*((float)x/TFT_WIDTH) - 1, 2*((float)y/TFT_HEIGHT) - 1}; // Convert screen pixels to (-1, 1) in x and y direction

    float alpha, beta, gamma;
    baycentric_coords(triag, v, alpha, beta, gamma); // Calculate baycentril coordinate of v for the given triangle.

    int time_1 = micros()-start_time;

    // If inside triangle
    if (alpha >= 0 && beta >= 0 && gamma >= 0) {
    //if (abs(alpha) < 0.01 || abs(beta) < 0.01 || abs(gamma) < 0.01) {

        Matrix<3> pos = mult(alpha, triag[0].pos) + mult(beta, triag[1].pos) + mult(gamma, triag[2].pos);
        Matrix<3> pos_s = mult(alpha, triag[0].pos_s) + mult(beta, triag[1].pos_s) + mult(gamma, triag[2].pos_s);

        int time_2 = micros()-start_time;

        // https://en.wikipedia.org/wiki/Z-buffering#Fixed-point_representation
        #ifdef Z_BUFF_16
            uint16_t z = (65536 - 1) * (f/(f-n) + (1/pos_s(2))*((-f*n)/(f-n))); 
        #else
            uint8_t z = (256 - 1) * (f/(f-n) + (1/pos_s(2))*((-f*n)/(f-n))); 
        #endif

        int time_3 = micros()-start_time;

        // If this new fragment is closer
        if (z < z_buffer[y-y_offset][x-x_offset]) {
            z_buffer[y-y_offset][x-x_offset] = z;

            Matrix<3> frag_color = mult(alpha, triag[0].color) + mult(beta, triag[1].color) + mult(gamma, triag[2].color);

            // Matrix<3> lm = norm(trans_light_pos - pos); // Light direction
            // Matrix<3> n = norm(mult(alpha, triag[0].norm) + mult(beta, triag[1].norm) + mult(gamma, triag[2].norm)); // Surface normal
            // Matrix<3> rm = norm(mult(2*(~lm*n)(0), n) - lm); // Reflected light
            // Matrix<3> v = norm(-pos); // viewer
            // Serial << "pos" << pos << "\n";
            // Serial << "light_pos" << trans_light_pos << "\n";
            // Serial << "lm" << lm << "\n";
            // Serial << "n" << n << "\n";
            // Serial << "rm" << rm << "\n";
            // Serial << "v" << v << "\n";

            int time_4 = micros()-start_time;

            // Matrix<3> ambient = mult(Ka, frag_color);
            // Matrix<3> diffuse = mult(Kd*max((~lm*n)(0), (float)0), frag_color);
            //Matrix<3> specular = mult(Ks*pow(max((~rm*v)(0), (float)0), SHININESS), light_color);
            Matrix<3> color = frag_color;// + specular;
        

            int time_5 = micros()-start_time;
            //Serial << color << "\n";

            //tft.draw_pixel(x, y, tft.matrix_to_565(color));

            // Mark pixel as shouldnt erase
            //should_erase[y][x/8] |= (0x01 << x%8);

            // if (frame_buffer[y][x] != color_332) {
            //     pixel_changed[y][x/8] |= (0x01 << x%8);
            // }

            #ifdef RGB565
                next_frame_buffer[y-y_offset][x-x_offset] = tft.matrix_to_565(color);
            #else
                next_frame_buffer[y-y_offset][x-x_offset] = tft.matrix_to_332(color);
            #endif

            int time_6 = micros()-start_time;

            //Serial << time_1 << " " << time_2-time_1 << " " << time_3-time_2 << " " << time_4-time_3 << " " << time_5-time_4 << " " << time_6-time_5 << "\n";
        }
    }
}

void Renderer::render() {
    int start_time = micros();

    int time_1 = micros()-start_time; // 964

    trans_light_pos = from_homo(view_matrix * light_pos); // Light position isnt transformed into screen cooridnate space. Kept in view space

    view_model_matrix = view_matrix * model_matrix;
    view_model_normal_matrix =~((view_model_matrix).Inverse());

    for (size_t i = 0; i < num_vertices; i++) {
        // Serial << "i" << (int)i << "\n";
        Matrix<3> color = {1,1,1};
        Vertex vertex = vertex_shader(dragon::vertices[i], dragon::normals[i/3], color);
        // Serial << "vertex pos" << vertex.pos << "\n";
        // Serial << "vertex norm" << vertex.norm << "\n";
        shaded_vertices[i] = vertex;
    }

    int time_2 = micros()-start_time; // 148
    
    int time_3 = micros()-start_time; // 6

    // For each block
    for (size_t block_x = 0; block_x < TFT_WIDTH/BLOCK_WIDTH; block_x++) {
        for (size_t block_y = 0; block_y < TFT_HEIGHT/BLOCK_HEIGHT; block_y++) {
    // for (size_t block_x = 0; block_x <= 0; block_x++) {
    //     for (size_t block_y = 0; block_y <= 0; block_y++) {
    
            // Clear depth and frame buffer
            for (size_t x = 0; x < BLOCK_WIDTH; x++) {
                for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
                    next_frame_buffer[y][x] = 0x00;
                    z_buffer[y][x] = INFINITY;
                    //frame_buffer[y][x] = 0x00;
                    //pixel_changed[y][x/8] = 0x00;
                    //should_erase[y][x/8] = 0xFF;
                }
            }

            uint32_t block_start_x = block_x*BLOCK_WIDTH;
            uint32_t block_start_y = block_y*BLOCK_HEIGHT;
            uint32_t block_end_x = min(block_start_x+BLOCK_WIDTH, (uint32_t)TFT_WIDTH);
            uint32_t block_end_y = min(block_start_y+BLOCK_HEIGHT, (uint32_t)TFT_HEIGHT);

            // For each triangle
            uint32_t min_x, min_y, max_x, max_y;
            for (size_t i = 0; i < num_indices; i+=3) {
                Vertex triag[3];
                triag[0] = shaded_vertices[i];//[dragon::indices[i]];
                triag[1] = shaded_vertices[i+1];//[dragon::indices[i+1]];
                triag[2] = shaded_vertices[i+2];//[dragon::indices[i+2]];

                if (is_backface(triag)) {
                    continue;
                }

                // Find bounding box and loop over those pixels
                triag_bbox(triag, min_x, min_y, max_x, max_y);

                // Clip to block bounds
                min_x = max(min_x, block_start_x);
                min_y = max(min_y, block_start_y);
                max_x = min(max_x, block_end_x-1);
                max_y = min(max_y, block_end_y-1);

                for (size_t x = min_x; x <= max_x; x++) {
                    for (size_t y = min_y; y <= max_y; y++) {
                        fragment_shader(x, y, block_start_x, block_start_y, triag);
                    }
                }
            }
            
            int time_4 = micros()-start_time; // 10151
            

            // for (size_t x = 0; x < TFT_WIDTH; x++) {
            //     for (size_t y = 0; y < TFT_HEIGHT; y++) {
            //         // If the pixel has not already been changed and is not black, we want to tell the program to make it black
            //         // if (!((pixel_changed[y][x/8] >> x%8) & 0x01) && (frame_buffer[y][x] != 0x00)) {
            //         //     pixel_changed[y][x/8] |= (0x01 << x%8);
            //         //     frame_buffer[y][x] = 0x00;
            //         // }
            //     }
            // }


            
            // for (size_t x = 0; x < TFT_WIDTH; x++) {
            //     for (size_t y = 0; y < TFT_HEIGHT; y++) {
                    
            //         // If pixel has changed since last frame, update it
            //         if ((pixel_changed[y][x/8] >> x%8) & 0x01) {
            //             tft.draw_pixel(x, y, frame_buffer[y][x]);
            //         }
                    
            //         //old_z_buffer[y][x] = z_buffer[y][x];
            //     }
            // }



            // Write differences between current frame buffer and next frame buffer
            // for (size_t x = block_start_x; x < block_end_x; x++) {
            //     for (size_t y = block_start_y; y < block_end_y; y++) {
                    
            //         //tft.draw_pixel(x, y, cur_frame_buffer[y/RES_REDUX][x/RES_REDUX]);
            //         //tft.draw_pixel(x, y, (z_buffer[y-block_start_y][x-block_start_x]>>10<<5));

            //         // // draw changes
            //         // if (cur_frame_buffer[y][x] != next_frame_buffer[y-block_start_y][x-block_start_x]) {    
            //         //     #ifdef RGB565
            //         //         //tft.draw_pixel(x, y, next_frame_buffer[y-block_start_y][x-block_start_x]);
            //         //     #else
            //         //         tft.draw_pixel(x, y, tft.RGB332_to_565(next_frame_buffer[y-block_start_y][x-block_start_x]));
            //         //     #endif
            //         // }
            //         // else {
            //         //     //tft.draw_pixel(x, y, 0xffff);
            //         // }

            //         cur_frame_buffer[y][x] = next_frame_buffer[y-block_start_y][x-block_start_x];
            //     }
            // }
            
            int time_5 = micros()-start_time; // 10522


            tft.draw_block(next_frame_buffer, block_start_x, block_start_y);

            // int res_redux_squared = RES_REDUX*RES_REDUX;
            // uint8_t r, g, b;
            // for (size_t x = 0; x < TFT_WIDTH/RES_REDUX; x++) {
            //     for (size_t y = 0; y < TFT_HEIGHT/RES_REDUX; y++) {
            //             uint16_t avg_r = 0;
            //             uint16_t avg_g = 0;
            //             uint16_t avg_b = 0;

            //             for (size_t i = 0; i < RES_REDUX; i++) {
            //                 for (size_t j = 0; j < RES_REDUX; j++) {
            //                     tft.RGB565_to_rgb(next_frame_buffer[RES_REDUX*y + j][RES_REDUX*x + i], r, g, b);
            //                     avg_r += r;
            //                     avg_g += g;
            //                     avg_b += b;
            //                 }
            //             }
                        
            //             avg_r /= res_redux_squared;
            //             avg_g /= res_redux_squared;
            //             avg_b /= res_redux_squared;
                        
            //             cur_frame_buffer[y][x] = tft.rgb_to_565(avg_r, avg_g, avg_b);

            //             //cur_frame_buffer[y][x] = next_frame_buffer[RES_REDUX*y][RES_REDUX*x];
            //             // random += 7;
            //     }
            // }
        }
    }
    
    int time_6 = micros()-start_time; // 10286
    
    Serial << "fps: " << (float)1e6/(micros()-start_time) << "\n";
    Serial << time_1 << " " << time_2-time_1 << " " << time_3-time_2 << " " << time_6-time_3 << "\n";
}

bool Renderer::is_backface(Vertex *triag) {
    Matrix<3> n = cross((triag[1].pos_s - triag[0].pos_s), (triag[2].pos_s - triag[0].pos_s)); // Normal of face

    return (~triag[0].pos_s * n)(0) < 0; 
}

float Renderer::mag(Matrix<4> &m) {
    return sqrtf(m(0)*m(0) + m(1)*m(1) + m(2)*m(2)) / m(3);
}

Matrix<4> Renderer::norm_homo(Matrix<4> &m) {
    Matrix<4> res = {m(0)/m(3), m(1)/m(3), m(2)/m(3), 1};
    return res;
}

Matrix<3> Renderer::from_homo(Matrix<4> const &m) {
    Matrix<3> res = {m(0)/m(3), m(1)/m(3), m(2)/m(3)};
    return res;
}

void Renderer::triag_bbox(Vertex *triag, uint32_t &min_x, uint32_t &min_y, uint32_t &max_x, uint32_t &max_y) {
    min_x = INFINITY; 
    min_y = INFINITY;
    max_x = -INFINITY; 
    max_y = -INFINITY;

    for (size_t i = 0; i < 3; i++) {
        // Convert from (-1, +1) to pixel coordinates
        float x = (triag[i].pos_s(0)/2 + (float)0.5) * TFT_WIDTH;
        float y = (triag[i].pos_s(1)/2 + (float)0.5) * TFT_HEIGHT;

        if (x < min_x) min_x = x;
        if (x > max_x) max_x = x;
        if (y < min_y) min_y = y;
        if (y > max_y) max_y = y;
    }
}
float Renderer::det(Matrix<2> u, Matrix<2> v) {
    return u(0)*v(1) - u(1)*v(0);
}

// Slide 107 Graphics lecture notes
void Renderer::baycentric_coords(Vertex *triag, Matrix<2> v, float &alpha, float &beta, float &gamma) {
    Matrix<2> a = {triag[0].pos_s(0), triag[0].pos_s(1)};
    Matrix<2> b = {triag[1].pos_s(0), triag[1].pos_s(1)};
    Matrix<2> c = {triag[2].pos_s(0), triag[2].pos_s(1)};

    static Vertex old_triag[3];
    static float x_delta_alpha, x_delta_beta, x_delta_gamma, y_delta_alpha, y_delta_beta, y_delta_gamma;
    static float alpha_0, beta_0, gamma_0;
    static bool deltas_computed;

    // If same as the last triangle
    if (triag[0] == old_triag[0] && triag[1] == old_triag[1] && triag[2] == old_triag[2]) {
        if (!deltas_computed) {
            Matrix<2> x_delta = {1,0};
            Matrix<2> y_delta = {0,1};
            Matrix<2> origin = {0,0};

            alpha_0 = ilf(c,b,origin) / ilf(c,b,a);
            beta_0  = ilf(a,c,origin) / ilf(a,c,b);
            gamma_0 = 1-alpha_0-beta_0;

            float alpha_x1 = ilf(c,b,x_delta) / ilf(c,b,a);
            float beta_x1  = ilf(a,c,x_delta) / ilf(a,c,b);
            float gamma_x1 = 1-alpha_x1-beta_x1;

            x_delta_alpha = alpha_x1 - alpha_0;
            x_delta_beta  = beta_x1  - beta_0;
            x_delta_gamma = gamma_x1 - gamma_0;

            float alpha_y1 = ilf(c,b,y_delta) / ilf(c,b,a);
            float beta_y1  = ilf(a,c,y_delta) / ilf(a,c,b);
            float gamma_y1 = 1-alpha_y1-beta_y1;

            y_delta_alpha = alpha_y1 - alpha_0;
            y_delta_beta  = beta_y1  - beta_0;
            y_delta_gamma = gamma_y1 - gamma_0;

            deltas_computed = true;
        }

        float x_delta = v(0);
        float y_delta = v(1);

        alpha = alpha_0 + x_delta*x_delta_alpha + y_delta*y_delta_alpha;
        beta  = beta_0  + x_delta*x_delta_beta  + y_delta*y_delta_beta;
        gamma = gamma_0 + x_delta*x_delta_gamma + y_delta*y_delta_gamma;
    }
    // If doesnt hit cache
    else {
        alpha = ilf(c,b,v) / ilf(c,b,a);
        beta  = ilf(a,c,v) / ilf(a,c,b);
        gamma = 1-alpha-beta;

        for (size_t i = 0; i < 3; i++) {
            old_triag[i] = triag[i];
        }
        
        deltas_computed = false;
    }
}

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @param v (-1, 1) vector to location on screen
 * @return float 
 */
float Renderer::ilf(Matrix<2> a, Matrix<2> b, Matrix<2> v) {
    float x = v(0);
    float y = v(1);

    return (a(1)-b(1))*x + (b(0)-a(0))*y + a(0)*b(1) - b(0)*a(1);
}

Matrix<3> Renderer::mult(float s, Matrix<3> const &m) {
    Matrix<3> res = {s*m(0), s*m(1), s*m(2)};

    return res;
}

Matrix<3> Renderer::norm(Matrix<3> const &m) {
    float mag = sqrtf(m(0)*m(0) + m(1)*m(1) + m(2)*m(2));
    Matrix<3> res = mult(1/mag, m);

    return res;
}

Matrix<3> Renderer::cross(Matrix<3> const &a, Matrix<3> const &b) {
    Matrix<3> res = {
        a(1)*b(2) - a(2)*b(1),
        -a(0)*b(2) + a(2)*b(0),
        a(0)*b(1) - a(1)*b(0)
    };

    return res;
}