#include <Arduino.h>
#include <SPI.h>
#include "TFT_Driver.h"
#include "TTGO_T_Display.h"
#include "ST7789.h"
#include <BasicLinearAlgebra.h>

SPIClass spi = SPIClass(VSPI);

uint8_t addr_col = 0;
uint8_t addr_row = 0;

bool locked = true;
bool inTransaction = false;
bool currently_is_command = false; 

void TFT_Driver::init() {
    spi.begin(TFT_SCLK, TFT_MISO, TFT_MOSI, -1);

    pinMode(TFT_CS, OUTPUT);
    digitalWrite(TFT_CS, HIGH); // Chip select high (inactive)

    pinMode(TFT_DC, OUTPUT);
    digitalWrite(TFT_DC, HIGH); // Data/Command high = data mode]

    pinMode(TFT_RST, OUTPUT);
    digitalWrite(TFT_RST, HIGH);
    delay(5);
    digitalWrite(TFT_RST, LOW);
    delay(20);
    digitalWrite(TFT_RST, HIGH);
    delay(150);
    write_command(TFT_SWRST);
    delay(150);

    begin_write();

    write_command(ST7789_SLPOUT);   // Sleep out
    delay(120);

    write_command(ST7789_NORON);    // Normal display mode on

    //------------------------------display and color format setting--------------------------------//
    write_command(ST7789_MADCTL);
    //write_data(0b11110000); 
    write_data(TFT_MAD_COLOR_ORDER); // Make bottom right (0,0)

    // JLX240 display datasheet
    write_command(0xB6);
    write_data(0x0A);
    write_data(0x82);

    write_command(ST7789_COLMOD);
    write_data(0x55);
    delay(10);

    //--------------------------------ST7789V Frame rate setting----------------------------------//
    write_command(ST7789_PORCTRL);
    write_data(0x0c);
    write_data(0x0c);
    write_data(0x00);
    write_data(0x33);
    write_data(0x33);

    write_command(ST7789_GCTRL);      // Voltages: VGH / VGL
    write_data(0x35);

    //---------------------------------ST7789V Power setting--------------------------------------//
    write_command(ST7789_VCOMS);
    write_data(0x28);		// JLX240 display datasheet

    write_command(ST7789_LCMCTRL);
    write_data(0x0C);

    write_command(ST7789_VDVVRHEN);
    write_data(0x01);
    write_data(0xFF);

    write_command(ST7789_VRHS);       // voltage VRHS
    write_data(0x10);

    write_command(ST7789_VDVSET);
    write_data(0x20);

    write_command(ST7789_FRCTR2);
    write_data(0x0f);

    write_command(ST7789_PWCTRL1);
    write_data(0xa4);
    write_data(0xa1);

    //--------------------------------ST7789V gamma setting---------------------------------------//
    write_command(ST7789_PVGAMCTRL);
    write_data(0xd0);
    write_data(0x00);
    write_data(0x02);
    write_data(0x07);
    write_data(0x0a);
    write_data(0x28);
    write_data(0x32);
    write_data(0x44);
    write_data(0x42);
    write_data(0x06);
    write_data(0x0e);
    write_data(0x12);
    write_data(0x14);
    write_data(0x17);

    write_command(ST7789_NVGAMCTRL);
    write_data(0xd0);
    write_data(0x00);
    write_data(0x02);
    write_data(0x07);
    write_data(0x0a);
    write_data(0x28);
    write_data(0x31);
    write_data(0x54);
    write_data(0x47);
    write_data(0x0e);
    write_data(0x1c);
    write_data(0x17);
    write_data(0x1b);
    write_data(0x1e);

    write_command(ST7789_INVON);

    write_command(ST7789_CASET);    // Column address set
    write_data(0x00);
    write_data(0x00);
    write_data(0x00);
    write_data(0xE5);    // 239

    write_command(ST7789_RASET);    // Row address set
    write_data(0x00);
    write_data(0x00);
    write_data(0x01);
    write_data(0x3F);    // 319

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    end_write();
    delay(120);
    begin_write();

    write_command(ST7789_DISPON);    //Display on
    delay(120);

    // Turn on the back-light LED
    pinMode(TFT_BL, OUTPUT);
    digitalWrite(TFT_BL, HIGH);

    end_write();
}

void TFT_Driver::begin_write() {
    if (locked) {
        locked = false;
        spi.beginTransaction(SPISettings(SPI_FREQUENCY, MSBFIRST, TFT_SPI_MODE));
        digitalWrite(TFT_CS, LOW); // Chip select low
    }
}
    
void TFT_Driver::end_write() {
    if(!inTransaction) {
        if (!locked) {
            locked = true;
            digitalWrite(TFT_CS, HIGH); // Chip select high
            spi.endTransaction();
        }
    }
}

void TFT_Driver::spi_write_8(uint8_t data, uint8_t is_command) {
    if (!is_command && currently_is_command) {
        DC_D;
        currently_is_command = false;
    }
    else if (is_command && !currently_is_command) {
        DC_C;
        currently_is_command = true;
    }

    spi.transfer(data);
}

void TFT_Driver::spi_write_16(uint16_t data, uint8_t is_command) {
    if (!is_command && currently_is_command) {
        DC_D;
        currently_is_command = false;
    }
    else if (is_command && !currently_is_command) {
        DC_C;
        currently_is_command = true;
    }
    
    spi.transfer16(data);
}

/**
 * @brief Send 8 bit command to display
 * 
 * @param data
 */
void TFT_Driver::write_command(uint8_t data) {
    //begin_write();
    //DC_C;
    spi_write_8(data, 1);
    //spi.transfer(data);
    //DC_D;
    //end_write();
}

void TFT_Driver::write_data(uint8_t data) {
    //begin_write();
    //DC_D;
    spi_write_8(data, 0);
    //spi.transfer(data);
    //digitalWrite(TFT_CS, LOW);; // "Allow more hold time for low VDI rail"?
    //end_write();
}

/**
 * @brief Draw single pixel
 * 
 * @param x 
 * @param y 
 * @param color 16 bit 565RGB color
 */
void TFT_Driver::draw_pixel(int32_t x, int32_t y, uint16_t color) {

    int start_time = micros();

    // Range checking
    if ((x < 0) || (y < 0) ||(x >= TFT_WIDTH) || (y >= TFT_HEIGHT)) return;

    x += COL_START;
    y += ROW_START;

    begin_write();

    // No need to send x if it has not changed (speeds things up)
    if (addr_col != x) {
        spi_write_8(TFT_CASET, 1);
        spi_write_16(x, 0); spi_write_16(x, 0); // Start and end x is the same
        addr_col = x;
    }

    // No need to send y if it has not changed (speeds things up)
    if (addr_row != y) {
        spi_write_8(TFT_RASET, 1);
        spi_write_16(y, 0); spi_write_16(y, 0); // Start and end y is the same
        addr_row = y;
    }

    // Write color info
    spi_write_8(TFT_RAMWR, 1);
    spi_write_16(color, 0);

    end_write();
}

void TFT_Driver::draw_buffer(uint16_t frame_buffer[TFT_HEIGHT][TFT_WIDTH]) {

    begin_write();

    set_window(0, 0, TFT_WIDTH-1, TFT_HEIGHT-1);

    spi_write_8(TFT_RAMWR, 1);
    for (size_t y = 0; y < TFT_HEIGHT; y++) {
        for (size_t x = 0; x < TFT_WIDTH; x++) {
            spi_write_16(frame_buffer[y][x], 0);   
        }
    }

    end_write();
}

void TFT_Driver::draw_block(uint16_t block_buffer[BLOCK_HEIGHT][BLOCK_WIDTH], uint32_t x_offset, uint32_t y_offset) {
    begin_write();

    set_window(x_offset, y_offset, x_offset + BLOCK_WIDTH-1, y_offset + BLOCK_HEIGHT-1);

    spi_write_8(TFT_RAMWR, 1);
    for (size_t y = 0; y < BLOCK_HEIGHT; y++) {
        for (size_t x = 0; x < BLOCK_WIDTH; x++) {
            spi_write_16(block_buffer[y][x], 0);   
        }
    }

    end_write();
}

void TFT_Driver::set_window(int32_t x0, int32_t y0, int32_t x1, int16_t y1) {
    //begin_tft_write(); // Must be called before setWindow

    addr_row = 0xFFFF;
    addr_col = 0xFFFF;

    x0 += COL_START;
    y0 += ROW_START;
    x1 += COL_START;
    y1 += ROW_START;

    // DC_C; 
    spi_write_8(TFT_CASET, 1);
    // DC_D; 
    spi_write_16(x0, 0); spi_write_16(x1, 0);
    // DC_C; 
    spi_write_8(TFT_RASET, 1);
    // DC_D; 
    spi_write_16(y0, 0); spi_write_16(y1, 0);

    //end_tft_write(); // Must be called after setWindow
}


void TFT_Driver::fill_rect(int32_t x, int32_t y, int32_t w, int32_t h, uint16_t color) {
  // Clipping
  if ((x >= TFT_WIDTH) || (y >= TFT_HEIGHT)) return;

  if (x < 0) { w += x; x = 0; }
  if (y < 0) { h += y; y = 0; }

  if ((x + w) > TFT_WIDTH)  w = TFT_WIDTH  - x;
  if ((y + h) > TFT_HEIGHT) h = TFT_HEIGHT - y;

  if ((w < 1) || (h < 1)) return;

  begin_write();

  set_window(x, y, x + w - 1, y + h - 1);

  write_color(color, w * h);

  end_write();
}

void TFT_Driver::write_color(uint16_t color, uint32_t len) {
    //DC_C; 
    spi_write_8(TFT_RAMWR, 1);
    
    //DC_D; 
    for (size_t i = 0; i < len; i++) {
        spi_write_16(color, 0);
    }
}

void TFT_Driver::fill_screen(uint16_t color) {
    fill_rect(0, 0, TFT_WIDTH, TFT_HEIGHT, color);
}

uint16_t TFT_Driver::rgb_to_565(uint8_t r, uint8_t g, uint8_t b) {
    uint8_t r5 = ( r * 249 + 1014 ) >> 11;
    uint8_t g6 = ( g * 253 +  505 ) >> 10;
    uint8_t b5 = ( b * 249 + 1014 ) >> 11;

    return (r5 << 11) | (g6 << 5) | b5;
}

uint16_t TFT_Driver::matrix_to_565(Matrix<3> const &m) {
    return rgb_to_565(min(m(0), (float)1)*255, min(m(1), (float)1)*255, min(m(2), (float)1)*255);
}

Matrix<3> TFT_Driver::RGB565_to_matrix(uint16_t color) {    
    uint8_t r5 = (color & 0b1111100000000000) >> 11;
    uint8_t g6 = (color & 0b0000011111100000) >> 5;
    uint8_t b5 = color & 0b0000000000011111;

    uint8_t r8 = ( r5 * 527 + 23 ) >> 6;
    uint8_t g8 = ( g6 * 259 + 33 ) >> 6;
    uint8_t b8 = ( b5 * 527 + 23 ) >> 6;

    Matrix<3> m = {r8/255.0, g8/255.0, b8/255.0};

    return m;
}

uint8_t TFT_Driver::rgb_to_332(uint8_t r, uint8_t g, uint8_t b) {
    uint8_t r3 = r >> 5;
    uint8_t g3 = g >> 5;
    uint8_t b2 = b >> 6;

    return (r3 << 5) | (g3 << 2) | b2;
}

uint8_t TFT_Driver::matrix_to_332(Matrix<3> const &m) {
    return rgb_to_332(min(m(0), (float)1)*255, min(m(1), (float)1)*255, min(m(2), (float)1)*255);
}

void TFT_Driver::RGB565_to_rgb(uint16_t color, uint8_t &r, uint8_t &g, uint8_t &b) {
    uint8_t r5 = (color & 0b1111100000000000) >> 11;
    uint8_t g6 = (color & 0b0000011111100000) >> 5;
    uint8_t b5 = color & 0b0000000000011111;

    r = ( r5 * 527 + 23 ) >> 6;
    g = ( g6 * 259 + 33 ) >> 6;
    b = ( b5 * 527 + 23 ) >> 6;
}

float TFT_Driver::color_delta(uint16_t a, uint16_t b) {

    Matrix<3> m_a = RGB565_to_matrix(a);
    Matrix<3> m_b = RGB565_to_matrix(b);

    Matrix<3> delta = m_b-m_a;

    float dist = sqrtf(powf(delta(0), 2) + powf(delta(1), 2) + powf(delta(2), 2));

    return dist/1.73205; // Normalize to 0-1
}

uint8_t TFT_Driver::RGB565_to_332(uint16_t rgb565) {
    uint8_t r3 = (rgb565 & 0b1111100000000000) >> 13;
    uint8_t g3 = (rgb565 & 0b0000011111100000) >> 8;
    uint8_t b2 = (rgb565 & 0b0000000000011111) >> 3;

    return (r3 << 5) | (g3 << 2) | b2; 
}

uint16_t TFT_Driver::RGB332_to_565(uint8_t RGB332) {
    uint8_t r5 = (RGB332 & 0b11100000) >> 2;
    uint8_t g6 = (RGB332 & 0b00011100) << 1;
    uint8_t b5 = (RGB332 & 0b00000011) << 3;

    return (r5 << 11) | (g6 << 5) | b5; 
}

uint8_t TFT_Driver::matrix_to_mono(Matrix<3> const &m) {
    return rgb_to_mono(min(m(0), (float)1)*255, min(m(1), (float)1)*255, min(m(2), (float)1)*255);
}

uint8_t TFT_Driver::rgb_to_mono(uint8_t r, uint8_t g, uint8_t b) {
    return (0.01328125 * r) + (0.0447125 * g) + (0.00450625 * b);
}

uint16_t TFT_Driver::mono_to_565(uint8_t mono) {
    return rgb_to_565(mono*16, mono*16, mono*16);
}