#version 330 core

layout (location = 0) in vec3 in_Pos;

out vec3 vertexPos;

void main()
{
    vertexPos = in_Pos;
    gl_Position = vec4(in_Pos.x, in_Pos.y, in_Pos.z, 1.0);
    gl_PointSize = 2.0;
}