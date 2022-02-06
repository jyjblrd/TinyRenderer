#include <SPI.h>
#include "TFT_Driver.h"
#include <BasicLinearAlgebra.h>
#include "Vertex.h"
#include "Renderer.h"

TFT_Driver tft = TFT_Driver();
Renderer renderer = Renderer(tft);

int start_time = millis();

void setup() {
  Serial.begin(2000000);

  tft.init();

  tft.fill_screen(0x00);

  // float n0[] = {0,0,1};
  // float n1[] = {0,-1,0};

  // float a0[] = {-1,-1,1};
  // float a1[] = {1,-1,1};
  // float a2[] = {-1,1,1};
  // float a3[] = {1,1,1};
  
  // float a4[] = {-1,-1,-1};
  // float a5[] = {1,-1,-1};

  // float a6[] = {-1,-2,1};
  // float a7[] = {-1,-2,-1};
  // float a8[] = {0,0,0};

  // float c0[] = {1,0,1};
  // float c1[] = {0,1,1};
  // float c2[] = {1,1,0};
  // Vertex v0 = Vertex(a0, n0, c0);
  // Vertex v1 = Vertex(a1, n0, c1);
  // Vertex v2 = Vertex(a2, n0, c2);

  // Vertex v3 = Vertex(a1, n0, c1);
  // Vertex v4 = Vertex(a3, n0, c0);
  // Vertex v5 = Vertex(a2, n0, c2);

  // Vertex v6 = Vertex(a0, n1, c0);
  // Vertex v7 = Vertex(a4, n1, c1);
  // Vertex v8 = Vertex(a5, n1, c2);

  // Vertex v9 = Vertex(a0, n1, c0);
  // Vertex v10 = Vertex(a5, n1, c2);
  // Vertex v11 = Vertex(a1, n1, c1);

  // Vertex v12 = Vertex(a6, n1, c0);
  // Vertex v13 = Vertex(a7, n1, c2);
  // Vertex v14 = Vertex(a8, n1, c1);

  // // Vertex vertices[] = {
  // //   v0, v1, v2, 
  // //   v3, v4, v5, 
  // //   v6, v7, v8, 
  // //   v9, v10, v11,
  // //   v12, v13, v14,
  // // };  
  // // int indices[] = {
  // //   0,1,2, 
  // //   3,4,5, 
  // //   6,7,8, 
  // //   9,10,11,
  // //   12,13,14,
  // // };
  
  // Vertex vertices[] = {
  //   v0, v1, v2, 
  //   v4, 
  //   v6, v7, v8, 
  //   v11,
  //   v12, v13, v14,
  // };
  // int indices[] = {
  //   0,1,2, 
  //   1,3,2, 
  //   4,5,6, 
  //   4,6,7,
  //   8,9,10,
  // };
  

  // renderer.set_vertices(dragon::vertices, dragon::normals, dragon::vertex_count);
  // renderer.set_indices(dragon::indices, dragon::vertex_count);

  renderer.calc_view_matrix();
}

void loop() {
  int elapsed_time = millis() - start_time;

  float theta = 2*PI * ((float)elapsed_time/10000); // Rotate once per 10 seconds

  Matrix<4,4> rotate_z = {
    cosf(theta),-sinf(theta),0,0,
    sinf(theta),cosf(theta),0,0,
    0,0,1,0,
    0,0,0,1
  };

  renderer.set_model_matrix(rotate_z);
  renderer.render();
}
