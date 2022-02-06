rendering methods:
- put eveyrthing in a frame buffer and then write that to the screen
- try keep track of which pixels need to change from last frame (didnt work, coudlnt keep track)
- keep previous frame in buffer but convert it to RGB332 so it uses less space (doesnt like black -> color transitions, also has uglier colors)
- make previous frame buffer a quater the size and compare to 2x2 pixel groups (looks nicer, but slightly(?)slower than 332 method)
- can try to compare color delta but too slow and doesnt work well
- render in blocks so that our next_frame_buffer is smaller and it can fit in ram. however this caused a huge ammount of screen tearing and looked like shit and wasnt that much faster
- Tried to use both cores by putting spi transfers into a buffer but the buffer had to be too but and used up too much memory so i scrapped it