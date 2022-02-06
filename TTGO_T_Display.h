#ifndef TTGO_T_Display_h
#define TTGO_T_Display_h

#define TFT_WIDTH  240
#define TFT_HEIGHT 135

#define BLOCK_WIDTH 80
#define BLOCK_HEIGHT 135

#define COL_START 40
#define ROW_START 53

#define TFT_MOSI            19
#define TFT_SCLK            18
#define TFT_MISO            -1
#define TFT_CS              5   // Chip select pin
#define TFT_DC              16  // Data/Command pin
#define TFT_RST             23

#define TFT_BL          4  // Display backlight control pin

#define SPI_FREQUENCY  60000000 // 6 MHz is the maximum SPI read speed for the ST7789V

#define DC_C digitalWrite(TFT_DC, LOW) // Command
#define DC_D digitalWrite(TFT_DC, HIGH) // Data

#define TFT_SPI_MODE SPI_MODE3

#endif