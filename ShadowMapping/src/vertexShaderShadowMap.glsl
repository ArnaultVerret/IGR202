#version 330 core            // minimal GL version support expected from the GPU

// input vertex data, different for all executions of this shader
layout(location=0) in vec3 vPositionModel;

// transformation bringing the scene into light's point of view
uniform mat4 depthMVP;
uniform mat4 model;

void main() {
  gl_Position = depthMVP * model * vec4(vPositionModel, 1);
}
