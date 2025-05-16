#version 330 core
in float vVal;
out vec4 FragColor;

void main() {
    //FragColor = vec4(vVal, vVal, vVal, vVal);
    
    FragColor = vec4(1,1,1,vVal);
}