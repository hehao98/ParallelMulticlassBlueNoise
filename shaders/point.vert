#version 330 core

layout (location = 0) in vec2 in_Pos;
layout (location = 1) in int  in_Type;

out vec2 vPos;
out vec4 vColor;

void main()
{
    vPos  = in_Pos;
    int vType = in_Type;
    if (vType == 0) {
        vColor = vec4(0.92, 0.23, 0.17, 1.0);
    } else if (vType == 1) {
        vColor = vec4(0.26, 0.58, 0.16, 1.0);
    } else if (vType == 2) {
        vColor = vec4(0.0, 0.0, 1.0, 1.0);
    } else if (vType == 3) {
        vColor = vec4(1.0, 1.0, 0.0, 1.0);
    } else if (vType == 4) {
        vColor = vec4(1.0, 0.0, 1.0, 1.0);
    } else if (vType == 5) {
        vColor = vec4(0.0, 1.0, 1.0, 1.0);
    } else {
        vColor = vec4(0.0, 0.0, 0.0, 1.0);
    }

    gl_Position = vec4(in_Pos.x, in_Pos.y, 0.0, 1.0);
    gl_PointSize = 4.0;
}