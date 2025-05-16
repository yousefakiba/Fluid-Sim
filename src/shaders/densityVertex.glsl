#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in float aVal;

out float vVal;

uniform mat4 uProjection;

void main() {
    vVal = aVal;
    gl_Position = uProjection * vec4(aPos, 0.0, 1.0);
    gl_PointSize = 20.0;
}