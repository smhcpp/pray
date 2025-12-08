#version 330
in vec2 fragTexCoord;
in vec4 fragColor;
out vec4 finalColor;
uniform sampler2D texture0;
uniform vec4 colDiffuse;
uniform vec2 player_pos;
uniform float radius;
uniform vec2 resolution;

void main()
{
    float mask = texture(texture0, fragTexCoord).r;
    vec2 pixelPos = fragTexCoord * resolution;
    vec2 corrected_player_pos = vec2(player_pos.x, resolution.y - player_pos.y);
    float dist = distance(pixelPos, corrected_player_pos);
    float minimum = 0.4;
    float falloff = clamp(1.0 - (dist / radius), minimum, 1.0);
    if (mask < minimum) mask = minimum;
    float visibility = mask * falloff;
    finalColor = vec4(0.0, 0.0, 0.0, 1.0 - visibility);
}
