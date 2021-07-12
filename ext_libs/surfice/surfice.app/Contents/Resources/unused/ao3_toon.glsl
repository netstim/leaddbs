//vert
#version 330
layout(location = 0) in vec3 Vert;
layout(location = 3) in vec2 Coord;
out vec2 texCoord;
void main () {
  gl_Position = vec4 (Vert, 1.0);
  texCoord = Coord;
}
//frag
#version 330
uniform sampler2D tex1, tex2, depth_texture1, depth_texture2;
uniform float blend1, alpha1, fracAO;
uniform vec2 texture_size;
smooth in vec2 texCoord;
out vec4 color;
vec2 x1 = vec2(1.0 / texture_size.x, 0.0);
vec2 y1 = vec2(0.0, 1.0 / texture_size.y);
float depth(in vec2 coord) {
	return texture(depth_texture1, coord).x;
}
float getDiff(vec2 xy) {

	float mz =  (depth(xy-x1-y1)+depth(xy-y1)+depth(xy+x1-y1)+depth(xy-x1)+depth(xy+x1)+depth(xy-x1+y1)+depth(xy+y1)+depth(xy+x1+y1));
	float dz = abs(depth(xy) - 0.125 * mz);
	dz = dz/mz * 8000.0;
	return smoothstep(0.0, 1.0, dz);
}
float getAO(void) {
	vec2 xy = texCoord;
	float dz = 0.125 *  (getDiff(xy-x1-y1)+getDiff(xy-y1)+getDiff(xy+x1-y1)+getDiff(xy-x1)+getDiff(xy+x1)+getDiff(xy-x1+y1)+getDiff(xy+y1)+getDiff(xy+x1+y1));
	return dz  * fracAO;
}
void main(void) {
  vec4 t1 = texture(tex1, texCoord);
  if (t1.a == 0.0) discard;
  vec4 t2 = texture(tex2, texCoord);
  if (fracAO > 0.0)
    t1.rgb = clamp(t1.rgb-getAO(), 0.0, 1.0);
  t1.rgb = mix(t2.rgb,t1.rgb, alpha1);
  float depth = 1.0 - (3.0 * (texture(depth_texture2, texCoord).x - texture(depth_texture1, texCoord).x));
  depth = clamp(depth, 0.0, 1.0);
  color = mix(t1, t2, blend1 * depth);
}