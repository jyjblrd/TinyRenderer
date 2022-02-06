#ifndef TFT_Driver_h
#define TFT_Driver_h

#include <Arduino.h>
#include <BasicLinearAlgebra.h>
#include "TTGO_T_Display.h"
using namespace BLA;

class TFT_Driver {
    public: 

        void init();

        void spi_write_8(uint8_t data, uint8_t is_command);
        void spi_write_16(uint16_t data, uint8_t is_command);

        void begin_write();
        void end_write();

        void write_command(uint8_t data);
        void write_data(uint8_t data);

        void draw_pixel(int32_t x, int32_t y, uint16_t color);
        void draw_buffer(uint16_t frame_buffer[TFT_HEIGHT][TFT_WIDTH]);
        void draw_block(uint16_t block_buffer[BLOCK_HEIGHT][BLOCK_WIDTH], uint32_t x_offset, uint32_t y_offset);

        void set_window(int32_t x0, int32_t y0, int32_t x1, int16_t y1);
        void fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint16_t color);
        void write_color(uint16_t color, uint32_t len);
        void fill_screen(uint16_t color);

        uint16_t rgb_to_565(uint8_t r, uint8_t g, uint8_t b);
        uint16_t matrix_to_565(Matrix<3> const &m);
        Matrix<3> RGB565_to_matrix(uint16_t color);
        uint8_t rgb_to_332(uint8_t r, uint8_t g, uint8_t b);
        uint8_t matrix_to_332(Matrix<3> const &m);
        float color_delta(uint16_t a, uint16_t b);
        uint8_t RGB565_to_332(uint16_t RGB565);
        uint16_t RGB332_to_565(uint8_t RGB332);
        void RGB565_to_rgb(uint16_t color, uint8_t &r, uint8_t &g, uint8_t &b);

        uint8_t matrix_to_mono(Matrix<3> const &m);
        uint8_t rgb_to_mono(uint8_t r, uint8_t g, uint8_t b);
        uint16_t mono_to_565(uint8_t mono);

        void run_queue();
};

#endif
