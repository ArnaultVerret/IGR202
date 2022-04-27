#version 330 core            // minimal GL version support expected from the GPU

uniform mat3 normMat;

struct LightSource {
  vec3 position;
  vec3 color;
  float intensity;
  int isActive;
  mat4 depthMVP;
};

int numberOfLights = 3;
uniform LightSource lightSources[3];
uniform sampler2D depthTex0;
uniform sampler2D depthTex1;
uniform sampler2D depthTex2;

struct Material {
  vec3 albedo;
  sampler2D albedoTex; // albedo map binded for the object
  int hasAlbedoTex; // auto albedo or mapping;
  sampler2D normalTex; // normal map binded for the object
  int hasNormalMap; // auto normals or mapping
};

uniform Material material;

uniform vec3 camPos;

in vec3 fPositionModel;
in vec3 fPosition;
in vec3 fNormal;
in vec2 fTexCoord;

in vec4 lightsPos[3];

out vec4 colorOut; // shader output: the color response attached to this fragment

float pi = 3.1415927;

float ShadowCalculation(LightSource light, int i){
  vec4 lightViewPos = lightsPos[i];
  vec3 projCoords = lightViewPos.xyz / lightViewPos.w;
  projCoords = projCoords * 0.5 + 0.5;
  float closestDepth = 0;
  if (i == 0)
    closestDepth = texture(depthTex0, projCoords.xy).r;
  else if (i == 1)
    closestDepth = texture(depthTex1, projCoords.xy).r;
  else
    closestDepth = texture(depthTex2, projCoords.xy).r;
  float currentDepth = projCoords.z;

  float shadow = currentDepth - 0.005 > closestDepth ? 1.0 : 0.0;

  return shadow;
}

// TODO: shadows
void main() {
  vec3 n = normalize(fNormal); // default normals
  if (material.hasNormalMap == 1){
    n = texture(material.normalTex, fTexCoord).rgb*2-1; // setup normals and put it in [-1, 1] range
  }
  vec3 albedo = material.albedo;
  if (material.hasAlbedoTex == 1){
    albedo = texture(material.albedoTex, fTexCoord).rgb;
  }
  vec3 wo = normalize(camPos - fPosition); // unit vector pointing to the camera

  vec3 radiance = vec3(0, 0, 0);
  for(int i=0; i<numberOfLights; ++i) {
    LightSource a_light = lightSources[i];
    if(a_light.isActive == 1) { // consider active lights only
      vec3 wi = normalize(a_light.position - fPosition); // unit vector pointing to the light
      vec3 Li = a_light.color*a_light.intensity;
      float shadow = ShadowCalculation(a_light, i);
      radiance += Li*albedo*max(dot(n, wi), 0)*(1.0-shadow);
    }
  }
  colorOut = vec4(radiance, 1.0); // build an RGBA value from an RGB one
}
